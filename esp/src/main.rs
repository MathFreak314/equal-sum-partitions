use esp::{MitmInfo, PRIMES, RunOpts, RunResult, crt2, factorial, run, run_mitm, sys};
use serde_json::json;
use std::io::Write;

const USAGE: &str = "usage: esp N [K] [--mitm] [--primes 0,1,2] [--record PATH] [--csv PATH] \
[--max-states X] [--mem-gb G] [--work-dir DIR] [--shard-bits B] [--threads T] [-q]";

fn utc_now() -> String {
    // Seconds since the epoch, formatted as ISO-8601 UTC without a date crate.
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs() as i64;
    let (days, sod) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + (mo <= 2) as i64;
    format!("{y:04}-{mo:02}-{d:02}T{:02}:{:02}:{:02}+00:00", sod / 3600, sod / 60 % 60, sod % 60)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let (mut primes, mut record, mut csv, mut max_states) = (vec![0usize], None, None, None);
    let (mut shard_bits, mut threads, mut verbose) = (12u32, 0usize, true);
    let (mut mem_gb, mut work_dir): (Option<f64>, Option<String>) = (None, None);
    let mut mitm = false;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| panic!("{a} needs a value\n{USAGE}")).clone();
        match a.as_str() {
            "--primes" => primes = val().split(',').map(|s| s.parse().expect("prime index")).collect(),
            "--record" => record = Some(val()),
            "--csv" => csv = Some(val()),
            "--max-states" => max_states = Some(val().parse().expect("max-states")),
            "--shard-bits" => shard_bits = val().parse().expect("shard-bits"),
            "--threads" => threads = val().parse().expect("threads"),
            "--mem-gb" => mem_gb = Some(val().parse().expect("mem-gb")),
            "--work-dir" => work_dir = Some(val()),
            "-q" | "--quiet" => verbose = false,
            "--mitm" => mitm = true,
            "-h" | "--help" => { println!("{USAGE}"); return; }
            _ => pos.push(a.parse::<u32>().unwrap_or_else(|_| panic!("bad argument {a}\n{USAGE}"))),
        }
    }
    assert!(!pos.is_empty() && pos.len() <= 2, "{USAGE}");
    assert!(primes.iter().all(|&i| i < PRIMES.len()), "prime indices are 0..{}", PRIMES.len() - 1);
    let n = pos[0] as usize;
    let k = pos.get(1).copied().unwrap_or(pos[0]);
    if threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    }
    let big_n = n as u32 * k;
    let avail = sys::avail_ram_bytes().unwrap_or(4 << 30);
    let mem_budget = mem_gb.map_or((avail as f64 * 0.75) as u64, |g| (g * (1u64 << 30) as f64) as u64);
    let work_dir = work_dir.map_or_else(std::env::temp_dir, std::path::PathBuf::from);
    println!("A({n},{k}): numbers 1..{big_n}, block sum {}, threads {}", k as f64 * (big_n as f64 + 1.0) / 2.0,
             rayon::current_num_threads());
    println!("memory budget {:.2} GB (available {:.2} GB), work dir {}", mem_budget as f64 / (1u64 << 30) as f64,
             avail as f64 / (1u64 << 30) as f64, work_dir.display());

    sys::keep_awake();
    let started = utc_now();
    let mut runs: Vec<(u64, RunResult)> = Vec::new();
    let mut infos: Vec<MitmInfo> = Vec::new();
    for &pi in &primes {
        let p = PRIMES[pi];
        println!("prime p{pi} = {p}{}", if mitm { " (meet in the middle)" } else { "" });
        let opts = RunOpts { fingerprint: record.is_some(), max_states, shard_bits, verbose, mem_budget,
                             work_dir: work_dir.clone(), ..RunOpts::new(n, k, p) };
        let r = if mitm {
            let (r, info) = run_mitm(&opts);
            println!("  mid states = {}  upper states = {}  dead mid states = {}  join passes = {}  join = {:.2}s",
                     info.mid_states, info.upper_states, info.dead_mid_states, info.join_passes, info.join_seconds);
            infos.push(info);
            r
        } else {
            run(&opts)
        };
        println!("  labeled mod p{pi} = {}{}  peak states = {}  time = {:.2}s  peak RSS = {:.1} MB  peak disk = {:.2} GB",
                 r.labeled_mod_p, if r.aborted { " (ABORTED)" } else { "" }, r.peak_states, r.total_seconds,
                 sys::peak_rss_bytes().unwrap_or(0) as f64 / 1048576.0, r.peak_disk_bytes as f64 / 1e9);
        runs.push((p, r));
    }
    let finished = utc_now();
    let aborted = runs.iter().any(|r| r.1.aborted);

    // Reconstruct the exact labeled count from the first two primes; check against the rest.
    let mut labeled: Option<u128> = None;
    let mut result: Option<u128> = None;
    if runs.len() >= 2 && !aborted {
        let x = crt2(runs[0].1.labeled_mod_p, runs[0].0, runs[1].1.labeled_mod_p, runs[1].0);
        let checks_ok = runs[2..].iter().all(|(p, r)| (x % *p as u128) as u64 == r.labeled_mod_p);
        let f = factorial(n);
        println!("labeled (CRT p{},p{}) = {x}", primes[0], primes[1]);
        println!("extra prime checks: {}", if runs.len() == 2 { "none".into() } else if checks_ok { "OK".to_string() } else { "MISMATCH".to_string() });
        println!("divisible by {n}! = {f}: {}", x % f == 0);
        if checks_ok && x % f == 0 {
            labeled = Some(x);
            result = Some(x / f);
            println!("result = {}", x / f);
        }
    }

    let first = &runs[0].1;
    if let Some(path) = csv {
        let mut f = std::fs::File::create(path).unwrap();
        writeln!(f, "m,states,transitions,seconds").unwrap();
        for l in &first.layers {
            writeln!(f, "{},{},{},{:.6}", l.m, l.states, l.transitions, l.seconds).unwrap();
        }
    }
    if let Some(path) = record {
        let layers: Vec<_> = first.layers.iter().enumerate().map(|(i, l)| {
            let col = |which: usize| -> Vec<String> {
                runs.iter().map(|(_, r)| r.layers.get(i).and_then(|x| x.fp)
                    .map_or("null".into(), |fp| if which == 0 { fp.0 } else { fp.1 }.to_string())).collect()
            };
            json!({
                "m": l.m, "states": l.states, "transitions": l.transitions,
                "seconds": runs.iter().map(|(_, r)| r.layers.get(i).map_or(0.0, |x| x.seconds)).collect::<Vec<_>>(),
                "key_bits": l.key_bits, "passes": l.passes,
                "input_on_disk": l.input_on_disk, "output_on_disk": l.output_on_disk,
                "count_sum_mod_p": col(0), "hash_sum_mod_p": col(1),
            })
        }).collect();
        let rec = json!({
            "schema": 1,
            "implementation": "rust",
            "code": {"crate": "esp", "version": env!("CARGO_PKG_VERSION"), "git_commit": env!("ESP_GIT_COMMIT"),
                     "git_dirty": env!("ESP_GIT_DIRTY")},
            "runtime": {"rustc": env!("ESP_RUSTC"), "threads": rayon::current_num_threads(),
                        "shard_bits": shard_bits, "mem_budget_bytes": mem_budget},
            "machine": {"os": std::env::consts::OS, "cpu": sys::cpu_brand(),
                        "logical_cpus": std::thread::available_parallelism().map(|x| x.get()).unwrap_or(0),
                        "ram_bytes": sys::total_ram_bytes()},
            "params": {"n_blocks": n, "block_size": k, "N": big_n, "S": k * (big_n + 1) / 2,
                       "max_states": max_states, "algorithm": if mitm { "mitm" } else { "forward" }},
            "mitm": infos.iter().map(|i| json!({
                "mid_states": i.mid_states, "upper_states": i.upper_states,
                "dead_mid_states": i.dead_mid_states, "join_passes": i.join_passes,
                "join_seconds": i.join_seconds})).collect::<Vec<_>>(),
            "primes": runs.iter().map(|(p, _)| p.to_string()).collect::<Vec<_>>(),
            "started_utc": started,
            "finished_utc": finished,
            "aborted": aborted,
            "result": result.map(|x| x.to_string()),
            "labeled": labeled.map(|x| x.to_string()),
            "labeled_mod_p": runs.iter().map(|(_, r)| r.labeled_mod_p.to_string()).collect::<Vec<_>>(),
            "peak_states": first.peak_states,
            "total_seconds": runs.iter().map(|(_, r)| r.total_seconds).collect::<Vec<_>>(),
            "peak_rss_bytes": sys::peak_rss_bytes(),
            "disk_bytes": runs.iter().map(|(_, r)| r.peak_disk_bytes).max(),
            "layers": layers,
        });
        std::fs::write(path, serde_json::to_string_pretty(&rec).unwrap() + "\n").unwrap();
    }
}
