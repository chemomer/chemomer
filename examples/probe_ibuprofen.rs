fn main() {
    let molfile = std::env::args().nth(1).expect("molfile path");
    let text = std::fs::read_to_string(&molfile).expect("read");
    let mut working = chemomer::parse_molfile(&text).expect("parse");
    chemomer::canonical::metals::materialize_hydrogens(&mut working);
    let sums = working.bond_order_sums();
    for (i, atom) in working.atoms.iter().enumerate() {
        let bs: Vec<u32> = working.bonds.iter().filter(|(a,b,_)| *a==i||*b==i).map(|(_,_,k)| *k).collect();
        println!("atom {} {} q={} h={} bonds={:?} sum={}", i + 1, atom.element, atom.charge, working.hydrogens()[i], bs, sums[i]);
    }
}
