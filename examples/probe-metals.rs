// quick probe via examples/
use chemomer::Engine;
fn main() {
    let cl = std::fs::read_to_string("/tmp/metalprobe/cl-.mol").unwrap();
    let nh2 = std::fs::read_to_string("/tmp/metalprobe/nh2-.mol").unwrap();
    for (name, text) in [("cl-", cl), ("nh2-", nh2)] {
        match Engine::from_molfile(&text) {
            Ok(e) => println!(
                "{name}: {} | {}",
                e.inchi_string()
                    .unwrap_or_else(|er| er.message().to_string()),
                e.aux_string().unwrap_or_default()
            ),
            Err(e) => println!("{name}: ERR {}", e.message()),
        }
    }
}
