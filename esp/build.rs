use std::process::Command;

fn run(cmd: &str, args: &[&str]) -> Option<String> {
    let out = Command::new(cmd).args(args).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn main() {
    let commit = run("git", &["rev-parse", "HEAD"]).unwrap_or_default();
    let dirty = run("git", &["status", "--porcelain", "--", "."]).map(|s| !s.is_empty());
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    println!("cargo:rustc-env=ESP_GIT_COMMIT={commit}");
    println!("cargo:rustc-env=ESP_GIT_DIRTY={}", dirty.map_or("unknown".into(), |d| d.to_string()));
    println!("cargo:rustc-env=ESP_RUSTC={}", run(&rustc, &["--version"]).unwrap_or_default());
    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=../.git/HEAD");
    println!("cargo:rerun-if-changed=../.git/index");
}
