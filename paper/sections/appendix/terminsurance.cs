calculation = 
{
  name = 'n-year term insurance',
  algorithm = { type = 'Runge Kutta 4', parameters = { stepsize = 0.01 } },
  equations = { 
      0 = { r_j = r, b_j = 0, mu_jk = { 1 = GM }, b_jk = { 1 = b01 } },
      1 = { },
  },
  range = { from = 50, to = 0 },
  boundaryvalues = { 0 = 0 , 1 = 0 },
  expressions = {
      interestrate = 0.05,
      bdeath = 1,
      n = 35,        
      age = 30,
      r(t) = interestrate,
      b01(t) = bdeath * indicator(t > 0) * indicator(t <= n),
      GM(t) = 0.0005 + 10 ^ (5.728 - 10 + 0.038*(age + t))
  }
}