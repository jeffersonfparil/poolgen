use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader, BufWriter};
use std::io::{Error, ErrorKind};
use std::sync::Arc;
use std::sync::Mutex;

pub fn parallel(fname: &String, n_threads: &u64, function: fn(u64, u64, ()) -> io::Result<String>, tuple_args: ()) -> io::Result<Arc<Mutex<Vec<String>>>> {
    // Find the positions whereto split the file into n_threads pieces
    let chunks = crate::io::parallel::find_file_splits(fname, n_threads);
    let n_digits = chunks[*n_threads as usize].to_string().len();
    println!("Chunks: {:?}", chunks);
    // Tuple arguments of read_chunks
    // Instantiate thread object for parallel execution
    let mut thread_objects = Vec::new();
    // Vector holding all returns from read_chunk()
    let mut thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..(*n_threads as usize) {
        let chunks_clone = chunks.clone();
        let tuple_args_clone = tuple_args.clone();
        // Determing start and end of this chunk
        let start = chunks_clone[i];
        let end = chunks_clone[i+1];
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let fname_out_per_thread = function(start, end, tuple_args_clone).unwrap();
            thread_ouputs_clone.lock().unwrap().push(fname_out_per_thread);
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }
    Ok(thread_ouputs)
}