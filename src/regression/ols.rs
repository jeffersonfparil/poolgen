use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use nalgebra::{self, DMatrix};
use std::sync::{Arc, Mutex};
use std::fs::{File, OpenOptions};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::sync::{Sync, sync_analyser_and_writer_single_thread};
use crate::io::sync::AlleleCountsOrFrequencies;
use crate::io::phen::{Phenotypes, load_phen};