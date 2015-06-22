extern crate num;

use num::bigint::{BigInt,ToBigInt};
use num::traits::{Zero,Signed,FromPrimitive};

use std::ops::{Neg,Mul,Sub,Add,Div};
use std::cmp::Ordering;
use std::vec::Vec;

struct Tableau {
    values: Vec<BigInt>,
    ncols: usize,
    nrows: usize,
    determinant: BigInt,
}

impl Tableau {
    fn new(ncols: usize, nrows: usize) -> Tableau {
        Tableau {
            values: vec![BigInt::zero(); ncols * nrows],
            ncols: ncols,
            nrows: nrows,
            determinant: BigInt::from_i32(-1).unwrap(),
        }
    }

    fn set(&mut self, row: usize, col: usize, value: BigInt) {
        self.values[row * self.ncols + col] = value;
    }

    fn entry(&self, row: usize, col: usize) -> &BigInt {
        &self.values[row * self.ncols + col]
    }

    fn set_det(&mut self, value: BigInt) {
        self.determinant = value;
    }

    fn pivot(&mut self, row: usize, col: usize) {

    	let pivelt = self.entry(row, col).clone(); /* pivelt anyhow later new determinant  */

    	if pivelt.is_zero() {
    		panic!("Trying to pivot on a zero");
    	}

        let cur_det = self.determinant.clone();
    	for i in 0..self.nrows {
    		if i != row {  // A[row][..] remains unchanged
                let nonzero = !self.entry(i, col).is_zero();
    			for j in 0..self.ncols {
    				if j != col {  // A[..][col] remains unchanged

    					//A[i,j] = (A[i,j] A[row,col] - A[i,col] A[row,j]) / det
    					let mut tmp1 = pivelt.abs().mul(self.entry(i, j));
    					if nonzero {
    						let tmp2 = self.entry(i, col).mul(self.entry(row, j));
    						if pivelt.is_negative() {
    							tmp1 = tmp1.add(tmp2);
    						} else {
    							tmp1 = tmp1.sub(tmp2);
    						}
    					}
    					self.set(i, j, tmp1.div(&cur_det));
    				}
    			}
    			if nonzero && !pivelt.is_negative() {
    				// row  i  has been dealt with, update  A[i][col] safely
                    let neg_entry = self.entry(i, col).neg();
    				self.set(i, col, neg_entry);
    			}
    		}
    	}

        self.set(row, col, cur_det);
    	if pivelt.is_negative() {
    		self.negate_row(row);
    	}

    	self.determinant = pivelt; //by construction always positive
    }

    fn negate_row(&mut self, row: usize) {
    	for j in 0..self.ncols {
    		let neg_entry = self.entry(row, j).neg();
    		self.set(row, j, neg_entry);
    	}
    }

    fn negate_col(&mut self, col: usize) {
    	for i in 0..self.nrows {
    		let neg_entry = self.entry(i, col).neg();
    		self.set(i, col, neg_entry);
    	}
    }

    // sign of  A[a,testcol] / A[a,col] - A[b,testcol] / A[b,col]
    // (assumes only positive entries of col are considered)
    fn ratio_test(&self, rowa: usize, rowb: usize, cola: usize, colb: usize) -> Ordering {
    	let a = self.entry(rowa, colb).mul(self.entry(rowb, cola));
    	let b = self.entry(rowb, colb).mul(self.entry(rowa, cola));
    	a.cmp(&b)
    }
}

#[test]
fn set_and_get_works() {

    let mut a = Tableau::new(5, 2);
    a.set(0, 0, BigInt::from_i32(2).unwrap());
    a.set(0, 1, BigInt::from_i32(2).unwrap());
    a.set(0, 2, BigInt::from_i32(1).unwrap());
    a.set(0, 3, BigInt::from_i32(-1).unwrap());
    a.set(1, 0, BigInt::from_i32(1).unwrap());
    a.set(1, 1, BigInt::from_i32(1).unwrap());
    a.set(1, 2, BigInt::from_i32(3).unwrap());
    a.set(1, 3, BigInt::from_i32(-1).unwrap());

    assert_eq!(&BigInt::from_i32(2).unwrap(), a.entry(0, 0));
    assert_eq!(&BigInt::from_i32(2).unwrap(), a.entry(0, 1));
    assert_eq!(&BigInt::from_i32(1).unwrap(), a.entry(0, 2));
    assert_eq!(&BigInt::from_i32(-1).unwrap(), a.entry(0, 3));
    assert_eq!(&BigInt::from_i32(1).unwrap(), a.entry(1, 0));
    assert_eq!(&BigInt::from_i32(1).unwrap(), a.entry(1, 1));
    assert_eq!(&BigInt::from_i32(3).unwrap(), a.entry(1, 2));
    assert_eq!(&BigInt::from_i32(-1).unwrap(), a.entry(1, 3));
}

#[test]
fn pivoting_works() {

    let n = 2;
	let mut a = Tableau::new(n, n+2);
	for i in 0..n {
		for j in 0..n+2 {
			let value = BigInt::from_usize((i + 1) + j*10).unwrap();
			a.set(i, j, value);
		}
	}

    assert_eq!(&BigInt::from_i32(1).unwrap(), a.entry(0, 0));
    assert_eq!(&BigInt::from_i32(11).unwrap(), a.entry(0, 1));
    assert_eq!(&BigInt::from_i32(2).unwrap(), a.entry(1, 0));
    assert_eq!(&BigInt::from_i32(12).unwrap(), a.entry(1, 1));

	a.pivot(0, 0);

    assert_eq!(&BigInt::from_i32(-1).unwrap(), a.entry(0, 0));
    assert_eq!(&BigInt::from_i32(11).unwrap(), a.entry(0, 1));
    assert_eq!(&BigInt::from_i32(-2).unwrap(), a.entry(1, 0));
    assert_eq!(&BigInt::from_i32(10).unwrap(), a.entry(1, 1));
}

#[test]
fn negating_col_works() {

	let n = 3;
	let mut a = Tableau::new(n, n+2);
	for i in 0..n {
		for j in 0..n+2 {
			let value = BigInt::from_usize(i + j*10).unwrap();
			a.set(i, j, value);
		}
	}

	a.negate_col(1);

    assert_eq!(&BigInt::from_i32(20).unwrap(), a.entry(0, 2));
    assert_eq!(&BigInt::from_i32(-10).unwrap(), a.entry(0, 1));
    assert_eq!(&BigInt::from_i32(-11).unwrap(), a.entry(1, 1));
    assert_eq!(&BigInt::from_i32(-12).unwrap(), a.entry(2, 1));
}

#[test]
fn positive_values_ratio_test_works() {

	let n = 2;
	let mut a = Tableau::new(n, n+2);
	for i in 0..n {
		for j in 0..n+2 {
			let value = BigInt::from_usize((i + 1) + j*10).unwrap();
			a.set(i, j, value);
		}
	}

    assert_eq!(Ordering::Greater, a.ratio_test(0, 1, 0, 1));
    assert_eq!(Ordering::Less, a.ratio_test(1, 0, 0, 1));
}
