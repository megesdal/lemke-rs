use super::tableau::Tableau;
use super::tableau_vars::{TableauVariable,TableauVariables};

use std::cmp::Ordering;
use num::traits::Signed;

#[cfg(test)] use num::bigint::BigInt;
#[cfg(test)] use num::traits::FromPrimitive;
#[cfg(test)] use time;

/*
 * minVar
 * ===========================================================
 * @return the leaving variable in  VARS, given by lexmin row,
 * when  enter  in VARS is entering variable
 * only positive entries of entering column tested
 * boolean  *z0leave  indicates back that  z0  can leave the
 * basis, but the lex-minratio test is performed fully,
 * so the returned value might not be the index of  z0
 */
pub fn lexminratio(tableau: &Tableau, vars: &TableauVariables, enter: &TableauVariable) -> (TableauVariable, bool) {

	let mut z0leave = false;
	let mut leave_candidate_rows = Vec::new();

	if vars.is_basic(&enter) {
		panic!("Variable {} is already in basis. Must be cobasic to enter.", enter.to_string());
	}

	let enter_col = vars.to_col(enter);

	// start with  leavecand = { i | A[i][col] > 0 }
	for i in 0..vars.n {
		if tableau.entry(i, enter_col).is_positive() {
			leave_candidate_rows.push(i);
		}
	}

	if leave_candidate_rows.len() == 0 {
		panic!("Ray termination when trying to enter {}", enter.to_string());  // convert to Result from panic...
	}

	/*else if (numcand == 1) {
		RecordStats(0, numcand);
		z0leave = IsLeavingRowZ0(leavecand[0]);
	}*/

	let z0_can_leave = process_candidates(tableau, vars, enter_col, &mut leave_candidate_rows);

	(vars.from_row(leave_candidate_rows[0]), z0leave)
}

/*
 * processCandidates
 * ================================================================
 * as long as there is more than one leaving candidate perform
 * a minimum ratio test for the columns of  j  in RHS, W(1),... W(n)
 * in the tableau.  That test has an easy known result if
 * the test column is basic or equal to the entering variable.
 */
fn process_candidates(tableau: &Tableau, vars: &TableauVariables, enter_col: usize, leave_candidate_rows: &mut Vec<usize>) -> bool {

	let z0_can_leave = process_rhs(tableau, vars, enter_col, leave_candidate_rows);
    let mut j = 1;
	while leave_candidate_rows.len() > 1 {
		//if j >= A.RHS() {                                             /* impossible, perturbed RHS should have full rank */
		//    throw new RuntimeException("lex-minratio test failed"); //TODO
		//}

		let wj = vars.w(j);
		if vars.is_basic(&wj) { /* testcol < 0: W(j) basic, Eliminate its row from leavecand */
			let removed = leave_candidate_rows.swap_remove(vars.to_row(&wj));
		} else { // not a basic testcolumn: perform minimum ratio tests
			let test_col = vars.to_col(&wj); /* since testcol is the  jth  unit column                    */
			if test_col != enter_col {   /* otherwise nothing will change */
				take_min_ratio_rows(tableau, enter_col, test_col, leave_candidate_rows);
			}
		}
        j += 1;
	}

	z0_can_leave
}

fn process_rhs(tableau: &Tableau, vars: &TableauVariables, enter_col: usize, leave_candidate_rows: &mut Vec<usize>) -> bool {

	take_min_ratio_rows(tableau, enter_col, vars.rhs_col(), leave_candidate_rows);

	let mut z0_can_leave = false;

	for i in 0..leave_candidate_rows.len() { // seek  z0  among the first-col leaving candidates
		if vars.from_row(leave_candidate_rows[i]).is_z0() {
            z0_can_leave = true;
			break;
		}
		/* alternative, to force z0 leaving the basis:
		 * return whichvar[leavecand[i]];
		 */
	}

	z0_can_leave
}

fn take_min_ratio_rows(tableau: &Tableau, enter_col: usize, test_col: usize, leave_candidate_rows: &mut Vec<usize>) {

	let mut num_min_candidates = 0;
	for i in 1..leave_candidate_rows.len() {  // investigate remaining candidates

		// sign of  A[l_0,t] / A[l_0,col] - A[l_i,t] / A[l_i,col]
		// note only positive entries of entering column considered
		match tableau.ratio_test(leave_candidate_rows[0], leave_candidate_rows[i], enter_col, test_col) {
            Ordering::Equal => {
                // new ratio is the same as before
    			num_min_candidates += 1;
    			leave_candidate_rows[num_min_candidates] = leave_candidate_rows[i];
            },
            Ordering::Greater => {
                // new smaller ratio detected
    			num_min_candidates = 0;
    			leave_candidate_rows[num_min_candidates] = leave_candidate_rows[i];
            },
            Ordering::Less => {}
        }
	}
    leave_candidate_rows.truncate(num_min_candidates+1);
}

/*
 * Test row candidate elimination using sign of:
 * A[l_0,t] / A[l_0,col] - A[l_i,t] / A[l_i,col]
 */
#[test]
fn take_min_ratio_elems_works() {

	let mut a = Tableau::new(2, 4);
	a.set(0, 0, BigInt::from_i32(2).unwrap());
	a.set(0, 1, BigInt::from_i32(2).unwrap());
	a.set(0, 2, BigInt::from_i32(1).unwrap());
    a.set(0, 3, BigInt::from_i32(-1).unwrap());
	a.set(1, 0, BigInt::from_i32(1).unwrap());
	a.set(1, 1, BigInt::from_i32(1).unwrap());
	a.set(1, 2, BigInt::from_i32(3).unwrap());
	a.set(1, 3, BigInt::from_i32(-1).unwrap());

	let mut candidates = vec![0, 1];

	let col1 = 1;
	let testcol1 = 2;
	let sgn1 = a.ratio_test(candidates[0], candidates[1], col1, testcol1);
	assert_eq!(Ordering::Less, sgn1);  // A[0,2] / A[0,1] - A[1,2] / A[1,1] should be 1/2 - 3/1 = -5/2

	take_min_ratio_rows(&a, col1, testcol1, &mut candidates);
	assert_eq!(1, candidates.len());
	assert_eq!(0, candidates[0]);

    candidates = vec![0, 1];
	let col2 = 2;
	let testcol2 = 1;
	let sgn2 = a.ratio_test(candidates[0], candidates[1], col2, testcol2);
	assert_eq!(Ordering::Greater, sgn2);  // A[0,1] / A[0,2] - A[1,1] / A[1,2] should be 2/1 - 1/3 = 5/3

	take_min_ratio_rows(&a, col2, testcol2, &mut candidates);
	assert_eq!(1, candidates.len());
	assert_eq!(1, candidates[0]);
}

#[test]
fn lexminvar_works() {

    let n = 2;
	let vars = TableauVariables::new(n);

    let mut a = Tableau::new(n, n+2);
	a.set(0, 0, BigInt::from_i32(2).unwrap());
	a.set(0, 1, BigInt::from_i32(2).unwrap());
	a.set(0, 2, BigInt::from_i32(1).unwrap());
    a.set(0, 3, BigInt::from_i32(-1).unwrap());
	a.set(1, 0, BigInt::from_i32(1).unwrap());
	a.set(1, 1, BigInt::from_i32(1).unwrap());
	a.set(1, 2, BigInt::from_i32(3).unwrap());
	a.set(1, 3, BigInt::from_i32(-1).unwrap());

	let (leave, z0_can_leave) = lexminratio(&a, &vars, &vars.z(0));
    assert_eq!(vars.w(2), leave);     // w2 = 4 is leaving
	assert_eq!(false, z0_can_leave);  // z0 can not leave

	let (leave, z0_can_leave) = lexminratio(&a, &vars, &vars.z(1));
    assert_eq!(vars.w(2), leave);     // w2 = 4 is leaving
	assert_eq!(false, z0_can_leave);  // z0 can not leave

	let (leave, z0_can_leave) = lexminratio(&a, &vars, &vars.z(2));
	assert_eq!(vars.w(1), leave);     // w1 = 3 is leaving
    assert_eq!(false, z0_can_leave);  // z0 can not leave

    // TODO: this plus ray termination...
	//_, _, err = lexminratio(lcp, lcp.w(1))
	//assert.NotNil(t, err, "Should be an error")
	//assert.Equal(t, "Variable w1 is already in basis. Must be cobasic to enter.", err.Error())

	//_, _, err = lexminratio(lcp, lcp.w(2))
	//assert.NotNil(t, err, "Should be an error")
	//assert.Equal(t, "Variable w2 is already in basis. Must be cobasic to enter.", err.Error())
}

#[test]
fn lexninvar_on_large_tableau_works() {

	let n = 1000;
    let vars = TableauVariables::new(n);

	let mut a = Tableau::new(n, n+2);

	for i in 0..a.nrows {
		for j in 0..a.ncols {
			if j == 0 {
				a.set(i, j, BigInt::from_i32(1).unwrap())
			} else {
                let v = {
                    let i = i as i64;
                    let j = j as i64;
                    (i-j+1)*(j*17-i*63)
                };
				a.set(i, j, BigInt::from_i64(v).unwrap());
			}
		}
	}

	let start = time::now();
	for i in 0..1000 {
		let (leave, z0_can_leave) = lexminratio(&a, &vars, &vars.z(0));
		assert_eq!(vars.w(1), leave);
		assert_eq!(false, z0_can_leave);
	}

	println!("1000 lexmin took: {}", time::now() - start);
}
