//! InChIKey — the truncated-SHA-256 key over the identifier, ported
//! from the reference's own generator (ikey_dll.c + ikey_base26.c;
//! MIT, source read at development time, never compiled or linked).
//!
//! Layout: the identifier splits after `InChI=1S/` into the MAJOR
//! layers (c, h, q — joined with their slashes), the optional PROTON
//! segment (/p±n), and the MINOR block (everything from the first
//! other layer). The major block hashes to four base-26 triplets plus
//! one dublet (14 chars); the minor block — concatenated with itself
//! once, the reference's quirk — hashes to two triplets plus one
//! dublet (8 chars); the blocks join with `-` and the tail carries
//! the calculation flag (S = standard), the subversion flag (A), a
//! dash, and the proton flag (N = none, else the ±1..±12 letter,
//! A beyond). The base-26 tables skip the E** and TAA..TTV triplets
//! and take the first 512 doublets.
//!
//! The 256-bit hash extensions (`XHash1`/`XHash2`, ikey_base26.c
//! `get_xtra_hash_{major,minor}_hex`) carry the digest bits the
//! truncated key did not consume: the major extension reads from the
//! first unconsumed bit (65) — 48 lowercase hex characters; the minor
//! extension from bit 37 — 56 characters. `-MergeHash` (ichimain.c)
//! renames the line to `InChIHash=` and dash-joins the requested
//! extensions behind the key.

use crate::error::{err, Result};

fn error_err(message: impl Into<String>) -> crate::error::Error {
    // err() returns Result; lift it to the error itself.
    err::<()>(message).unwrap_err()
}

/// The InChIKey for a standard identifier.
pub fn inchikey(inchi: &str) -> Result<String> {
    key_extensions(inchi).map(|(key, _, _)| key)
}

/// The key plus both 256-bit hash extensions: `(standard key, XHash1,
/// XHash2)`. The extensions ride the same SHA-256 digests the
/// truncated key reads — they are the un-consumed digest bits rendered
/// in lowercase hex (see [`xhash_major_hex`] / [`xhash_minor_hex`]).
pub fn key_extensions(inchi: &str) -> Result<(String, String, String)> {
    let blocks = split_blocks(inchi)?;
    let (digest_major, digest_minor) = digests(&blocks);
    let key = assemble_key(blocks.proton.as_deref(), &digest_major, &digest_minor);
    Ok((
        key,
        xhash_major_hex(&digest_major),
        xhash_minor_hex(&digest_minor),
    ))
}

/// The XHash1 extension of the FIRST key block (the reference's
/// `get_xtra_hash_major_hex`, ikey_base26.c): the major digest bits
/// the 65-bit truncated key did not consume — byte 8 with bit 64 (the
/// dublet's last bit) masked off `0xfe`, then bytes 9..31 — as 48
/// lowercase hex characters (191 bits).
pub fn xhash_major_hex(digest: &[u8; 32]) -> String {
    let mut out = format!("{:02x}", digest[8] & 0xfe);
    for byte in &digest[9..] {
        out.push_str(&format!("{byte:02x}"));
    }
    out
}

/// The XHash2 extension of the SECOND key block (the reference's
/// `get_xtra_hash_minor_hex`, ikey_base26.c): the minor digest bits
/// the 37-bit truncated key did not consume — byte 4 with bits 32..36
/// (the dublet's tail) masked off `0xe0`, then bytes 5..31 — as 56
/// lowercase hex characters (219 bits).
pub fn xhash_minor_hex(digest: &[u8; 32]) -> String {
    let mut out = format!("{:02x}", digest[4] & 0xe0);
    for byte in &digest[5..] {
        out.push_str(&format!("{byte:02x}"));
    }
    out
}

/// The MergeHash composition (ichimain.c): the key line renames to
/// `InChIHash=` and each requested extension dash-joins behind the
/// key — an absent extension leaves a bare dash in its slot.
pub fn merged_line(key: &str, xtra1: Option<&str>, xtra2: Option<&str>) -> String {
    format!(
        "InChIHash={}-{}-{}",
        key,
        xtra1.unwrap_or(""),
        xtra2.unwrap_or("")
    )
}

/// The identifier's split blocks: the MAJOR layers (formula, c, h, q
/// — joined with their slashes), the MINOR block (everything from the
/// first other layer, leading slash included — the reference's
/// `sminor` copy starts at the block's own slash), and the PROTON
/// segment (`/p±n`).
struct Blocks {
    major: String,
    minor: String,
    proton: Option<String>,
    /// The calculation flag character: S standard, N full (the
    /// FixedH /f layer).
    flag: char,
}

fn split_blocks(inchi: &str) -> Result<Blocks> {
    // Both identifier flavors hash through the same tables; the flag
    // character carries which: S standard, N full (the FixedH /f
    // layer, ichiparm.c "FixedH").
    let (body, flag) = if let Some(rest) = inchi.strip_prefix("InChI=1S/") {
        (rest, 'S')
    } else if let Some(rest) = inchi.strip_prefix("InChI=1/") {
        (rest, 'N')
    } else {
        return error_err("inchikey: not an InChI=1[S/] identifier");
    };
    let body = body.to_string();
    let segments: Vec<&str> = body.split('/').collect();
    let mut major: Vec<&str> = Vec::new();
    let mut proton: Option<&str> = None;
    let mut minor: Vec<&str> = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        if !minor.is_empty() {
            minor.push(segment);
            continue;
        }
        if index == 0 {
            // The formula segment leads the major block.
            major.push(segment);
            continue;
        }
        let letter = segment.chars().next().unwrap_or(' ');
        match letter {
            'c' | 'h' | 'q' => major.push(segment),
            'p' => proton = Some(segment),
            _ => minor.push(segment),
        }
    }
    if major.is_empty() {
        return err("inchikey: identifier has no major layers");
    }
    let smajor = major.join("/");
    // The reference's minor copy starts at the minor block's own
    // slash (ikey_dll.c copies str + j with j at the slash).
    let sminor = if minor.is_empty() {
        String::new()
    } else {
        format!("/{}", minor.join("/"))
    };
    Ok(Blocks {
        flag,
        major: smajor,
        minor: sminor,
        proton: proton.map(str::to_string),
    })
}

/// The truncated SHA-256 digests over both blocks. The minor block —
/// concatenated with itself once (the reference's quirk,
/// ikey_dll.c:410–415) — hashes to the minor digest.
fn digests(blocks: &Blocks) -> ([u8; 32], [u8; 32]) {
    let digest_major = sha256(blocks.major.as_bytes());
    let mut minor_bytes = blocks.minor.as_bytes().to_vec();
    if !minor_bytes.is_empty() && minor_bytes.len() < 255 {
        minor_bytes.extend_from_slice(&minor_bytes.clone());
    }
    let digest_minor = sha256(&minor_bytes);
    (digest_major, digest_minor)
}

fn assemble_key(proton: Option<&str>, digest_major: &[u8; 32], digest_minor: &[u8; 32]) -> String {
    let mut key = String::new();
    for chunk in 0..4 {
        key.push_str(&triplet(digest_major, chunk * 14));
    }
    key.push_str(&dublet(digest_major, 56));
    key.push('-');
    key.push_str(&triplet(digest_minor, 0));
    key.push_str(&triplet(digest_minor, 14));
    key.push_str(&dublet(digest_minor, 28));
    key.push(blocks.flag); // the calculation flag: S standard, N full
    key.push('A'); // the base subversion flag
    key.push('-');
    key.push_str(&proton_flag(proton));
    key
}

/// The proton flag: N for none, else the ±1..±12 letters, A beyond.
fn proton_flag(proton: Option<&str>) -> String {
    let Some(segment) = proton else {
        return "N".to_string();
    };
    // strtol-style over everything after the 'p' — the sign rides
    // along; unparseable tails read as 0.
    let value: i32 = segment[1..].trim().parse().unwrap_or(0);
    match value {
        0 => "N".to_string(),
        1..=12 => ["O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"][value as usize - 1]
            .to_string(),
        -12..=-1 => ["M", "L", "K", "J", "I", "H", "G", "F", "E", "D", "C", "B"]
            [(-value) as usize - 1]
            .to_string(),
        _ => "A".to_string(),
    }
}

/// The `index`-th base-26 triplet of the digest: 14 bits of the digest
/// bit stream (bytes in order, bits LSB-first within a byte) against
/// the triplet table that skips the E** and TAA..TTV entries. Shared
/// with the RInChI key shapes (rinchi_hashing.cpp uses the same
/// tables).
pub(crate) fn triplet(digest: &[u8; 32], index: usize) -> String {
    let mut value: usize = 0;
    for bit in 0..14 {
        let position = index + bit;
        let byte = digest[position / 8];
        if byte & (1 << (position % 8)) != 0 {
            value |= 1 << bit;
        }
    }
    let _ = value;
    table_entry(value)
}

/// The `index`-th base-26 dublet: 9 bits (512 entries). Shared with
/// the RInChI key shapes.
pub(crate) fn dublet(digest: &[u8; 32], index: usize) -> String {
    let mut value: usize = 0;
    for bit in 0..9 {
        let position = index + bit;
        let byte = digest[position / 8];
        if byte & (1 << (position % 8)) != 0 {
            value |= 1 << bit;
        }
    }
    let first = char::from_u32(u32::try_from(value / 26).expect("dublet first") + b'A' as u32)
        .expect("letter");
    let second = char::from_u32(u32::try_from(value % 26).expect("dublet second") + b'A' as u32)
        .expect("letter");
    format!("{first}{second}")
}

/// The v-th triplet of AAA..ZZZ in lexicographic order with the E**
/// block and the TAA..TTV block removed — 16384 entries, matching the
/// reference's tabulation (kept T entries start at TTW).
pub fn table_entry(value: usize) -> String {
    let mut skipped = 0;
    let mut candidate = [b'A'; 3];
    loop {
        let starts_e = candidate[0] == b'E';
        let in_tail_block = candidate[0] == b'T'
            && (candidate[1] < b'T' || (candidate[1] == b'T' && candidate[2] <= b'V'));
        if !starts_e && !in_tail_block {
            if skipped == value {
                return String::from_utf8_lossy(&candidate).into_owned();
            }
            skipped += 1;
        }
        // Lexicographic advance.
        for position in (0..3).rev() {
            if candidate[position] < b'Z' {
                candidate[position] += 1;
                break;
            }
            candidate[position] = b'A';
        }
    }
}

/// SHA-256 (FIPS 180-4) — dependency-free, the crate carries no
/// runtime dependencies. Shared with the RInChI key shapes.
pub(crate) fn sha256(message: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut state: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let mut padded = message.to_vec();
    let bit_length = (message.len() as u64) * 8;
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_length.to_be_bytes());

    for block in padded.chunks(64) {
        let mut w = [0u32; 64];
        for (index, word) in block.chunks(4).enumerate() {
            w[index] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for index in 16..64 {
            let s0 = w[index - 15].rotate_right(7)
                ^ w[index - 15].rotate_right(18)
                ^ (w[index - 15] >> 3);
            let s1 = w[index - 2].rotate_right(17)
                ^ w[index - 2].rotate_right(19)
                ^ (w[index - 2] >> 10);
            w[index] = w[index - 16]
                .wrapping_add(s0)
                .wrapping_add(w[index - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h) = (
            state[0], state[1], state[2], state[3], state[4], state[5], state[6], state[7],
        );
        for index in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let choose = (e & f) ^ ((!e) & g);
            let temp1 = h
                .wrapping_add(s1)
                .wrapping_add(choose)
                .wrapping_add(K[index])
                .wrapping_add(w[index]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let majority = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(majority);
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }
        state[0] = state[0].wrapping_add(a);
        state[1] = state[1].wrapping_add(b);
        state[2] = state[2].wrapping_add(c);
        state[3] = state[3].wrapping_add(d);
        state[4] = state[4].wrapping_add(e);
        state[5] = state[5].wrapping_add(f);
        state[6] = state[6].wrapping_add(g);
        state[7] = state[7].wrapping_add(h);
    }
    let mut digest = [0u8; 32];
    for (index, word) in state.iter().enumerate() {
        digest[index * 4..index * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_matches_the_published_vectors() {
        let empty = sha256(b"");
        assert_eq!(
            empty.iter().map(|b| format!("{b:02x}")).collect::<String>(),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        let abc = sha256(b"abc");
        assert_eq!(
            abc.iter().map(|b| format!("{b:02x}")).collect::<String>(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn the_triplet_table_skips_e_and_the_t_tail() {
        // v=0 -> AAA; the E** block is skipped (DZZ at 2703, FAA at
        // 2704); the T tail keeps only TTW..TZZ (TWA sits at 12224,
        // the first kept T entry).
        assert_eq!(table_entry(0), "AAA");
        assert_eq!(table_entry(1), "AAB");
        assert_eq!(table_entry(25), "AAZ");
        assert_eq!(table_entry(26), "ABA");
        assert_eq!(table_entry(2703), "DZZ");
        assert_eq!(table_entry(2704), "FAA");
        assert_eq!(table_entry(12224), "TWA");
        assert_eq!(table_entry(12172), "TUA");
        assert_eq!(table_entry(12515), "UHF");
    }

    #[test]
    fn water_keys_like_the_oracle() {
        // The recorded identity/water fixture. (The widely-published
        // pre-1.07 key differs in the hash blocks; the corpus records
        // what the 1.07.5 oracle produces.)
        let key = inchikey("InChI=1S/H2O/h1H2").expect("keys");
        assert_eq!(key, "XLYOFNOQVPJJNP-UHFFFAOYSA-N");
    }

    #[test]
    fn l_alanine_keys_like_the_oracle() {
        // The recorded identity/l-alanine fixture (proton flag N even
        // with a /p-free string; the /t /m /s layers ride the minor).
        let key =
            inchikey("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1").expect("keys");
        assert_eq!(key, "QNAYBMKLOCPYGJ-UWTATZPHSA-N");
    }

    #[test]
    fn the_extensions_match_the_oracle_for_water() {
        // The clean 1.07.5 oracle's `-Key -XHash1 -XHash2` output for
        // InChI=1S/H2O/h1H2: the major digest over "H2O/h1H2", the
        // minor digest over the empty block (sha256("")). The minor
        // extension's leading `80` is byte 4 of e3b0c44298.. masked
        // 0xe0 — the empty-minor shape, visible right in the recording.
        let (key, xhash1, xhash2) = key_extensions("InChI=1S/H2O/h1H2").expect("keys");
        assert_eq!(key, "XLYOFNOQVPJJNP-UHFFFAOYSA-N");
        assert_eq!(xhash1.len(), 48);
        assert_eq!(xhash1, "e014f429aafac730dc565eecf4ec4d2c4263e05c04ae2ef6");
        assert_eq!(xhash2.len(), 56);
        assert_eq!(
            xhash2,
            "80fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn the_merged_line_composes_like_the_oracle() {
        // ichimain.c: the line renames to InChIHash= and the requested
        // extensions dash-join; an absent extension leaves a bare dash.
        let line = merged_line(
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            Some("e014f429aafac730dc565eecf4ec4d2c4263e05c04ae2ef6"),
            Some("80fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        );
        assert_eq!(
            line,
            "InChIHash=XLYOFNOQVPJJNP-UHFFFAOYSA-N-\
             e014f429aafac730dc565eecf4ec4d2c4263e05c04ae2ef6-\
             80fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            merged_line("XLYOFNOQVPJJNP-UHFFFAOYSA-N", None, None),
            "InChIHash=XLYOFNOQVPJJNP-UHFFFAOYSA-N--"
        );
        assert_eq!(
            merged_line("XLYOFNOQVPJJNP-UHFFFAOYSA-N", Some("ab"), None),
            "InChIHash=XLYOFNOQVPJJNP-UHFFFAOYSA-N-ab-"
        );
    }

    #[test]
    fn the_extension_lengths_follow_the_reference_comments() {
        // ikey_base26.c: major extension 256-65 = 191 bits, minor
        // 256-37 = 219 bits — hex-rendered as 48 and 56 characters.
        let (_, xhash1, xhash2) =
            key_extensions("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1")
                .expect("keys");
        assert_eq!(xhash1.len(), 48);
        assert_eq!(xhash2.len(), 56);
        assert!(xhash1
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()));
        assert!(xhash2
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()));
    }
}
