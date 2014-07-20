calculation = 
{
  name = 'Disability annuity',
  algorithm = { type = 'Runge Kutta 4', parameters = { stepsize = 0.01 } },
  equations = { 
      0 = { r_j = r, b_j = 0, mu_jk = { 1 = GM01, 2 = GM02 }, b_jk = { } },
      1 = { r_j = r, b_j = b1, mu_jk = { 2 = GM12 }, b_jk = { } },
      2 = { },
  },
  range = { from = 50, to = 0 },
  boundaryvalues = { 0 = 0 , 1 = 0 , 2 = 0 },
  expressions = {
      interestrate = 0.05,
      bdisabled = 1,
      n = 35,        
      age = 30,
      r(t) = interestrate,
      b1(t) = bdisabled * indicator(t > 0) * indicator(t < n),
      GM01(t) = 0.0006 + 10 ^ (4.71609 - 10 + 0.06*(age + t)),
      GM02(t) = 0.0005 + 10 ^ (5.728 - 10 + 0.038*(age + t)),
      GM12(t) = GM02(t)
  }
}