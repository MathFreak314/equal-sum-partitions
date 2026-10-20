//! Machine info for run records: peak RSS, total RAM, CPU brand string.

#[cfg(windows)]
mod win {
    #[repr(C)]
    #[derive(Default)]
    pub struct ProcessMemoryCounters {
        pub cb: u32,
        pub page_fault_count: u32,
        pub peak_working_set_size: usize,
        pub working_set_size: usize,
        pub quota_peak_paged_pool_usage: usize,
        pub quota_paged_pool_usage: usize,
        pub quota_peak_non_paged_pool_usage: usize,
        pub quota_non_paged_pool_usage: usize,
        pub pagefile_usage: usize,
        pub peak_pagefile_usage: usize,
    }

    #[repr(C)]
    #[derive(Default)]
    pub struct MemoryStatusEx {
        pub length: u32,
        pub memory_load: u32,
        pub total_phys: u64,
        pub avail_phys: u64,
        pub total_page_file: u64,
        pub avail_page_file: u64,
        pub total_virtual: u64,
        pub avail_virtual: u64,
        pub avail_extended_virtual: u64,
    }

    #[link(name = "kernel32")]
    unsafe extern "system" {
        pub fn GetCurrentProcess() -> isize;
        pub fn K32GetProcessMemoryInfo(h: isize, p: *mut ProcessMemoryCounters, cb: u32) -> i32;
        pub fn GlobalMemoryStatusEx(p: *mut MemoryStatusEx) -> i32;
        pub fn SetThreadExecutionState(flags: u32) -> u32;
    }
}

/// Ask Windows not to sleep while this process runs (released automatically on exit).
pub fn keep_awake() {
    #[cfg(windows)]
    unsafe {
        const ES_CONTINUOUS: u32 = 0x8000_0000;
        const ES_SYSTEM_REQUIRED: u32 = 0x0000_0001;
        win::SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);
    }
}

pub fn peak_rss_bytes() -> Option<u64> {
    #[cfg(windows)]
    {
        let mut c = win::ProcessMemoryCounters::default();
        c.cb = size_of::<win::ProcessMemoryCounters>() as u32;
        let ok = unsafe { win::K32GetProcessMemoryInfo(win::GetCurrentProcess(), &mut c, c.cb) };
        (ok != 0).then_some(c.peak_working_set_size as u64)
    }
    #[cfg(not(windows))]
    None
}

/// (total, available) physical memory.
fn ram_bytes() -> Option<(u64, u64)> {
    #[cfg(windows)]
    {
        let mut m = win::MemoryStatusEx::default();
        m.length = size_of::<win::MemoryStatusEx>() as u32;
        let ok = unsafe { win::GlobalMemoryStatusEx(&mut m) };
        (ok != 0).then_some((m.total_phys, m.avail_phys))
    }
    #[cfg(not(windows))]
    None
}

pub fn total_ram_bytes() -> Option<u64> {
    ram_bytes().map(|r| r.0)
}

pub fn avail_ram_bytes() -> Option<u64> {
    ram_bytes().map(|r| r.1)
}

pub fn cpu_brand() -> String {
    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::x86_64::__cpuid;
        let mut bytes = Vec::with_capacity(48);
        for leaf in 0x8000_0002u32..=0x8000_0004 {
            let r = __cpuid(leaf);
            for w in [r.eax, r.ebx, r.ecx, r.edx] {
                bytes.extend_from_slice(&w.to_le_bytes());
            }
        }
        String::from_utf8_lossy(&bytes).trim_matches(|c: char| c == '\0' || c.is_whitespace()).to_string()
    }
    #[cfg(not(target_arch = "x86_64"))]
    "unknown".to_string()
}
