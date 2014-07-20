calculation = 
{
  name = 'Pure endowment',
  algorithm = { type = 'Runge Kutta 4', parameters = { stepsize = 0.01 } },
  equations = { 
      0 = { r_j = r, b_j = b0, mu_jk = { 1 = GM }, b_jk = { } },
      1 = { },
  },
  range = { from = 40, to = 0 },
  boundaryvalues = { 0 = 0 , 1 = 0 },
  expressions = {
      interestrate = 0.05,
      bpension = 1,
      pensiontime = 35,
      age = 30,
      r(t) = interestrate,
      b0(t) = bpension * delta(t - pensiontime),
      GM(t) = 0.0005 + 10 ^ (5.728 - 10 + 0.038*(age + t))
  }
}