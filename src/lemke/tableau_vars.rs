use num::bigint::BigInt;
use num::rational::{Ratio,BigRational};
use num::traits::{One,Zero};

use std::fmt::{Formatter,Debug,Error};
use std::ops::{Div,Rem,Mul};

use super::tableau::Tableau;

pub struct TableauVariable {
    value: usize,
    n: usize,
}

impl TableauVariable {

    pub fn is_z(&self) -> bool {
    	self.value <= self.n
    }

    pub fn is_w(&self) -> bool {
    	self.value > self.n
    }

    pub fn is_z0(&self) -> bool {
        self.value == 0
    }

    /*
     * complement of  v  in VARS, error if  v==Z(0).
     * this is  W(i) for Z(i)  and vice versa, i=1...n
     */
    pub fn complement(&self) -> TableauVariable {

    	if self.value == 0 {
    		panic!("Attempt to find complement of z0.");
    	}

    	if self.is_z() {
            TableauVariable{ value: self.value + self.n, n: self.n }
    	} else {
            TableauVariable{ value: self.value - self.n, n: self.n }
        }
    }
}

impl ToString for TableauVariable {

    fn to_string(&self) -> String {
        if self.is_z() {
    	    format!("z{}", self.value)
    	} else {
    	    format!("w{}", self.value - self.n)
        }
    }
}

impl PartialEq for TableauVariable {
    fn eq(&self, other: &TableauVariable) -> bool {
        self.value == other.value
    }
}

impl Debug for TableauVariable {
    fn fmt(&self, formatter: &mut Formatter) -> Result<(), Error> {
        formatter.write_str(&self.to_string())
    }
}

pub struct TableauVariables {
    /*  v in VARS, v cobasic:  TABCOL(v) is v's tableau col */
	/*  v  basic:  TABCOL(v) < 0,  TABCOL(v)+n   is v's row */
	/* VARS   = 0..2n = Z(0) .. Z(n) W(1) .. W(n)           */
	/* ROWCOL = 0..2n,  0 .. n-1: tabl rows (basic vars)    */
	/*                  n .. 2n:  tabl cols  0..n (cobasic) */
	vars2rowcol: Vec<usize>,
	rowcol2vars: Vec<usize>,
    pub n: usize,
}

impl TableauVariables {
    /*
     * init tableau variables:
     * Z(0)...Z(n)  nonbasic,  W(1)...W(n) basic
     * This is for setting up a complementary basis/cobasis
     */
    pub fn new(n: usize) -> TableauVariables {

        let mut vars = TableauVariables {
            n: n,
            vars2rowcol: vec![0; 2*n+1],
            rowcol2vars: vec![0; 2*n+1],
        };

        // w[1..n]
        for i in 1..n+1 {
            vars.rowcol2vars[i-1] = i + n;
            vars.vars2rowcol[i+n] = i - 1;
        }

        // z[0..n]
        for i in 0..n+1 {
    		vars.rowcol2vars[i+n] = i;
            vars.vars2rowcol[i] = n + i;
    	}

        vars
    }

    // TODO: check bounds with panic?
    pub fn z(&self, idx: usize) -> TableauVariable {
        TableauVariable{ value: idx, n: self.n }
    }

    pub fn w(&self, idx: usize) -> TableauVariable {
        TableauVariable{ value: idx + self.n, n: self.n }
    }

    pub fn from_row(&self, row: usize) -> TableauVariable {
        TableauVariable{ value: self.rowcol2vars[row], n: self.n }
    }

    pub fn from_col(&self, col: usize) -> TableauVariable {
    	TableauVariable{ value: self.rowcol2vars[col+self.n], n: self.n }
    }

    pub fn to_row(&self, var: &TableauVariable) -> usize {
        self.vars2rowcol[var.value]
    }

    pub fn to_col(&self, var: &TableauVariable) -> usize {
        self.vars2rowcol[var.value] - self.n
    }

    pub fn is_basic(&self, var: &TableauVariable) -> bool {
    	self.vars2rowcol[var.value] < self.n
    }

    pub fn rhs_col(&self) -> usize {
        self.n + 1
    }

    pub fn negate_rhs(&self, tableau: &mut Tableau) {
        tableau.negate_col(self.rhs_col())
    }

    fn swap(&mut self, enter: &TableauVariable, leave: &TableauVariable) -> (usize, usize) {

    	let leave_row = self.to_row(leave); // basic var is leaving
    	let enter_col = self.to_col(enter); // cobasic var is entering

    	self.vars2rowcol[leave.value] = enter_col + self.n;
    	self.rowcol2vars[enter_col+self.n] = leave.value;

    	self.vars2rowcol[enter.value] = leave_row;
    	self.rowcol2vars[leave_row] = enter.value;

    	(leave_row, enter_col)
    }

    /*
     * Pivot tableau on the element  A[row][col] which must be nonzero
     * afterwards tableau normalized with positive determinant
     * and updated tableau variables
     * @param leave (r) VAR defining row of pivot element
     * @param enter (s) VAR defining col of pivot element
     */
     pub fn pivot(&mut self, tableau: &mut Tableau, leave: &TableauVariable, enter: &TableauVariable) {

     	if !self.is_basic(leave) {
     		panic!("{} is not in the basis", leave.to_string());
     	}

     	if self.is_basic(enter) {
     		panic!("{} is already in the basis", enter.to_string());
     	}

     	let (row, col) = self.swap(enter, leave);  // update tableau variables
     	println!("pivoting ({},{})", row, col);
     	tableau.pivot(row, col)
     }

     /*
      * LCP result
      * current basic solution turned into  solz [0..n-1]
      * note that Z(1)..Z(n)  become indices  0..n-1
      */
     pub fn solution(&self, tableau: &Tableau, scale_factors: &Vec<BigInt>) -> Vec<BigRational> {

     	let mut z = Vec::with_capacity(self.n);
     	for i in 1..self.n + 1 {
     		z.push(self.result(tableau, scale_factors, &self.z(i)))
     	}
     	z
     }

     /*
      * Z(i):  scfa[i]*rhs[row] / (scfa[RHS]*det)
      * W(i):  rhs[row] / (scfa[RHS]*det)
      */
     fn result(&self, tableau: &Tableau, scale_factors: &Vec<BigInt>, var: &TableauVariable) -> BigRational {
     	if self.is_basic(var) {

            let one = BigInt::one();
            let row = self.to_row(var);
     		let scale_factor = if var.is_z() {
     			&scale_factors[row]
     		} else {
     			&one
     		};

            let col = self.rhs_col();
     		let numer = scale_factor.mul(tableau.entry(row, col));
     		let denom = (&tableau.determinant).mul(&scale_factors[col]);

     		Ratio::new(numer, denom)
     	} else {
             BigRational::zero()
         }
     }
}

#[test]
fn variable_assignments_work() {

	let vars = TableauVariables::new(4);

	for i in 0..vars.n+1 {
		let zi = vars.z(i);
		assert_eq!(i, zi.value);  // z(i) should equal i
		assert_eq!(i, vars.to_col(&zi));  // col(z(i)) should equal i
		assert_eq!(zi, vars.from_col(i));  // var(col(z(i))) = var(i) should equal z(i)
		assert_eq!(false, vars.is_basic(&zi));  // z(i) should NOT be basic
	}

	for i in 1..vars.n+1 {
		let wi = vars.w(i);
		assert_eq!(i+vars.n, wi.value);  // w(i) should equal i + n
		assert_eq!(i-1, vars.to_row(&wi));  // row(w(i)) should equal i - 1 (w is 1-indexed)
		assert_eq!(wi, vars.from_row(i-1));  // var(row(w(i - 1))) = var(i - 1) should equal w(i)
		assert_eq!(true, vars.is_basic(&wi));  // w(i) should be basic
	}
}

#[test]
fn complement_works() {

    let vars = TableauVariables::new(4);

	for i in 1..vars.n+1 {
		let comp_var = vars.z(i).complement();
		assert_eq!(vars.w(i), comp_var);  // w(i) should be complement of z(i)
	}

	//_, err := lcp.complement(lcp.z(0))
	//assert.NotNil(t, err, "Should not be able to get complement of z0")
}

#[test]
fn swap_works() {

	let mut vars = TableauVariables::new(4);

	let leave_var = vars.w(1);                    // first row, w(1) = 4
	assert_eq!(true, vars.is_basic(&leave_var));    // leaving var should be basic
	assert_eq!(0, vars.to_row(&leave_var));        // w(1) points to row 0

	let enter_var = vars.z(0);                    // first col, z(0) = 0
	assert_eq!(false, vars.is_basic(&enter_var));  // entering var should be cobasic
	assert_eq!(0, vars.to_col(&enter_var));        // z(0) points to col 0

	let (row, col) = vars.swap(&enter_var, &leave_var);

	assert_eq!(0, col);  // col should be original col(w1)
	assert_eq!(0, row);  // row should be original row(z0)

	assert_eq!(leave_var, vars.from_col(col));  // w1 should be var(col)
	assert_eq!(enter_var, vars.from_row(row));  // z0 should be var(row)
}
