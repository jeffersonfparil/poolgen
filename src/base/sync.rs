use std::io::{self, prelude::*, Error, ErrorKind, BufReader, BufWriter, SeekFrom};
use std::fs::{File, OpenOptions};
use std::str;
use nalgebra::DMatrix;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};
use crate::base::*;

// Parse a line of pileup into PileupLine struct
