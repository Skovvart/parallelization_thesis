// Simple C# implementation of Runge-Kutta 4 solvers for Thiele's
// differential equation, with Actulus calculation kernel examples.
// sestoft@itu.dk * 2012-01-11, 2012-02-01, 2013-01-29

// 2013-02-28: Added a hack to avoid the discontinuities at entire
// years.  This considerably improves convergence, so that RK4 at 100
// steps/year produces very accurate results.

// DO NOT DISTRIBUTE: Contains project-internal information.

// The primary purpose of this code is pedagogical; hence the
// increasingly complicated methods RK4_1 --> RK4_2 --> RK4_n.

// In the end, only the most general function is needed:
// double[][] RK4_n(dV, a, b, steps, Va) 
// where dV of type Action<double, double[], double[]> implements 
// an arbitrary number of derivatives.  This one is used to run all 
// the Actulus calculation specification examples.

// The weird encapsulation style (classes PureEndowment and so on) is
// motivated by the wish to stay relatively close to C, so no
// inheritance.  Also, while the function argument dV could be
// implemented in C using a function pointer, this is probably
// impossible or a bad idea on a GPU.  Instead the various definitions
// of the b_j and b_jk and mu_jk functions etc could simply be in
// separate source files, each with its own copy of RK4_n.

using System;
using System.Text;
using System.Diagnostics;

class CalculationSpecifications
{
    static readonly int steps = 100; // Per year

    static readonly double age = 30,
      interestrate = 0.05,
      bpension = 1,
      pensiontime = 35;

    static double indicator(bool b)
    {
        return b ? 1.0 : 0.0;
    }

    const double ln10 = 2.3025850929940459;

    // Gompertz-Makeham mortality intensities for Danish women
    static double GM(double t)
    {
        return 0.0005 + Math.Exp(ln10 * (5.728 - 10 + 0.038 * (age + t)));
    }

    static double r(double t)
    {
        return interestrate;    // Fixed interest rate
        // return rFsa(t);            // Finanstilsynet's rate curve
    }

    // The Danish FSA yield curve (Finanstilsynets rentekurve).
    // Data from 2011-11-16 
    static readonly double[]
      ts = new double[] { 
      0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 },
      rs = new double[] { 
      1.146677033, 1.146677033, 1.146677033, 1.340669678, 1.571952911, 1.803236144, 
      2.034519377, 2.26580261, 2.497085843, 2.584085843, 2.710085843, 2.805085843, 
      2.871485843, 2.937885843, 3.004285843, 3.070685843, 3.137085843, 3.136485843, 
      3.135885843, 3.135285843, 3.134685843, 3.134085843, 3.113185843, 3.092285843, 
      3.071385843, 3.050485843, 3.029585843, 3.008685843, 2.987785843, 2.966885843, 
      2.945985843, 2.925085843
    };

    // Get discount rate at time t by linear interpolation into (ts,rs); then
    // compute the instantaneous forward rate as described in 
    // https://wiki.actulus.dk/Documentation-CalculationPlatform-YieldCurves.ashx

    // This method uses binary search, which is needlessly slow because
    // the t values are monotonically decreasing.  Hence it would be
    // faster to keep the most recent index m into the ts/rs arrays and
    // decrement m only when t < ts[m].  It would also be easier to get
    // wrong, because it relies on many assumptions, so for now we stick
    // to the binary search.

    static double rFsa(double t)
    {
        // Requires ts non-empty and elements strictly increasing.
        int last = ts.Length - 1;
        if (t <= ts[0])
            return Math.Log(1 + rs[0] / 100);
        else if (t >= ts[last])
            return Math.Log(1 + rs[last] / 100);
        else
        {
            int a = 0, b = last;
            // Now a < b (bcs. ts must have more than 1 element) and ts[a] < t < ts[b]
            while (a + 1 < b)
            {
                // Now a < b and ts[a] <= t < ts[b]
                int i = (a + b) / 2;
                if (ts[i] <= t)
                    a = i;
                else // t < ts[i]
                    b = i;
            }
            // Now a+1>=b and ts[a] <= t < ts[b]; so a!=b and hence a+1 == b <= last
            int m = a;
            double tm = ts[m], tm1 = ts[m + 1];
            double rm = rs[m] / 100, rm1 = rs[m + 1] / 100;
            double Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm);
            return Math.Log(1 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1 + Rt);
        }
    }

    // Payment stream in state 0 (alive)
    static double b_0(double t)
    {
        return 0.0;
    }

    // Lump sum payments while in state 0 (alive) 
    static double bj_00(double t)
    {
        // This works only because t is known to be an integer
        return t == pensiontime ? bpension : 0.0;
    }

    // Lump sum payments while in state 1 (dead) 
    static double bj_11(double t)
    {
        // This works only because t is known to be an integer
        return 0.0;
    }

    // Transition intensity from state 0 (alive) to state 1 (dead)
    static double mu_01(double t)
    {
        return GM(t);
    }

    // Lump sum payment on transition from state 0 (alive) to state 1 (dead)
    static double bj_01(double t)
    {
        return 0.0;
    }

    // The dV0/dt function as used in Thiele:

    static double dV0(double t, double V0)
    {
        return r(t) * V0 - b_0(t) - mu_01(t) * (0 - V0 + bj_01(t));
    }

    static double dV0(double t, double V0, double V1)
    {
        return r(t) * V0 - b_0(t) - mu_01(t) * (0 - V0 + bj_01(t));
    }

    static double dV1(double t, double V0, double V1)
    {
        return 0.0;
    }

    static double[] dV(double t, double[] V)
    {
        double dV0 = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t));
        double dV1 = 0;
        return new double[] { dV0, dV1 };
    }

    // Fixed-step one-state Runge-Kutta solver RK4_1.
    // Solve from a>=b years in steps per year starting with V0(a)=Va.

    static double[] RK4_1(int a, int b, int steps, double Va)
    {
        double[] result = new double[a - b + 1];
        double h = -1.0 / steps;
        result[a - b] = Va;
        for (int y = a; y > b; y--)
        {
            double v = result[y - b];
            v += bj_00(y);
            double t = y;
            for (int s = 0; s < steps; s++)
            { 	// Integrate over [y,y+1]
                double k1 = h * dV0(t, v);
                double k2 = h * dV0(t + h / 2, v + k1 / 2);
                double k3 = h * dV0(t + h / 2, v + k2 / 2);
                double k4 = h * dV0(t + h, v + k3);
                v += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                t += h;
            }
            result[y - 1 - b] = v;
        }
        return result;
    }

    // Fixed-step two-state Runge-Kutta solver RK4_2.
    // Solve from a>=b years in steps per year starting with Vi(a)=Via.

    static double[][] RK4_2(int a, int b, int steps, double V0a, double V1a)
    {
        double[][] result = new double[a - b + 1][];
        for (int y = a; y >= b; y--)
            result[y - b] = new double[2];
        double h = -1.0 / steps;
        result[a - b][0] = V0a;
        result[a - b][1] = V1a;
        for (int y = a; y > b; y--)
        {
            double v0 = result[y - b][0], v1 = result[y - b][1];
            v0 += bj_00(y);
            v1 += bj_11(y);
            double t = y;
            for (int s = 0; s < steps; s++)
            { 	// Integrate over [y,y+1]
                double k01 = h * dV0(t, v0, v1);
                double k11 = h * dV1(t, v0, v1);
                double k02 = h * dV0(t + h / 2, v0 + k01 / 2, v1 + k11 / 2);
                double k12 = h * dV1(t + h / 2, v0 + k01 / 2, v1 + k11 / 2);
                double k03 = h * dV0(t + h / 2, v0 + k02 / 2, v1 + k12 / 2);
                double k13 = h * dV1(t + h / 2, v0 + k02 / 2, v1 + k12 / 2);
                double k04 = h * dV0(t + h, v0 + k03, v1 + k13);
                double k14 = h * dV1(t + h, v0 + k03, v1 + k13);
                v0 += (k01 + 2 * k02 + 2 * k03 + k04) / 6;
                v1 += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
                t += h;
            }
            result[y - 1 - b][0] = v0;
            result[y - 1 - b][1] = v1;
        }
        return result;
    }

    // Fixed-step n-state Runge-Kutta RK4_n solver.
    // Solve from a>=b years in steps per year starting with V0(a)=Va.
    // Note: It is very inefficient to create new k1, k2, k3, k4 arrays 
    // in every step.  In a C implementation all four could be allocated 
    // outside the loops, and sax and saxpy should write to those
    // preallocated arrays.  This would obscure the code somewhat.

    static double[][] RK4_n(int a, int b, int steps, double[] Va)
    {
        int n = Va.Length;
        double[][] result = new double[a - b + 1][];
        for (int y = a; y >= b; y--)
            result[y - b] = new double[n];
        double h = -1.0 / steps;
        result[a - b] = Va;
        for (int y = a; y > b; y--)
        {
            double[] v = result[y - b];
            double t = y;
            for (int s = 0; s < steps; s++)
            { 	// Integrate backwards over [y,y-1]
                double[] k1 = sax(h, dV(t, v));
                double[] k2 = sax(h, dV(t + h / 2, saxpy(0.5, k1, v)));
                double[] k3 = sax(h, dV(t + h / 2, saxpy(0.5, k2, v)));
                double[] k4 = sax(h, dV(t + h, saxpy(1, k3, v)));
                v = saxpy(1 / 6.0, k1, saxpy(2 / 6.0, k2, saxpy(2 / 6.0, k3, saxpy(1 / 6.0, k4, v))));
                t += h;
            }
            Array.Copy(v, result[y - 1 - b], v.Length);
        }
        return result;
    }

    // General n-state version, works for any set V of equations.  Slow
    // "declarative" but understandable version without memory reuse.

    static double[][] RK4_n(Func<double, double[], double[]> dV,
                            Func<double, double[]> bj_ii,
                            int a, int b, int steps, double[] Va)
    {
        int n = Va.Length;
        double[][] result = new double[a - b + 1][];
        for (int y = a; y >= b; y--)
            result[y - b] = new double[n];
        double h = -1.0 / steps;
        result[a - b] = Va;
        for (int y = a; y > b; y--)
        {
            double[] v = result[y - b];
            v = saxpy(1.0, v, bj_ii(y));
            double t = y;
            for (int s = 0; s < steps; s++)
            {     // Integrate backwards over [y,y-1]
                double[] k1 = sax(h, dV(t, v));
                double[] k2 = sax(h, dV(t + h / 2, saxpy(0.5, k1, v)));
                double[] k3 = sax(h, dV(t + h / 2, saxpy(0.5, k2, v)));
                double[] k4 = sax(h, dV(t + h, saxpy(1.0, k3, v)));
                v = saxpy(1 / 6.0, k1, saxpy(2 / 6.0, k2, saxpy(2 / 6.0, k3, saxpy(1 / 6.0, k4, v))));
                t += h;
            }
            Array.Copy(v, result[y - 1 - b], v.Length);
        }
        return result;
    }

    // sax = scalar a times x array, "declarative" version
    static double[] sax(double a, double[] x)
    {
        double[] res = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
            res[i] = a * x[i];
        return res;
    }

    // saxpy = scalar a times x array plus y array, "declarative" version
    static double[] saxpy(double a, double[] x, double[] y)
    {
        if (x.Length != y.Length)
            throw new Exception("saxpy: lengths of x and y differ");
        double[] res = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
            res[i] = a * x[i] + y[i];
        return res;
    }

    // General n-state version, works for any set V of equations.  Fast
    // "imperative" but obscure version that reuses intermediate arrays.

    static double[][] RK4_n(Action<double, double[], double[]> dV,
                Action<double, double[]> bj_ii,
                int a, int b, int steps, double[] Va)
    {
        int n = Va.Length;
        double[][] result = new double[a - b + 1][];
        for (int y = a; y >= b; y--)
            result[y - b] = new double[n];
        double h = -1.0 / steps;
        result[a - b] = Va;
        double[]
          k1 = new double[n],
          k2 = new double[n],
          k3 = new double[n],
          k4 = new double[n],
          tmp = new double[n];
        double[] v = new double[n];
        for (int y = a; y > b; y--)
        {
            bj_ii(y, v);
            v = saxpy(1.0, v, result[y - b]);
            for (int s = 0; s < steps; s++)
            { 	// Integrate backwards over [y,y-1]
                double t = y - s / (double)steps;
                // Hack: Fake limit from left
                dV(s == 0 ? t - 1e-14 : t, v, k1);
                sax(h, k1, k1);
                saxpy(0.5, k1, v, tmp);
                dV(t + h / 2, tmp, k2);
                sax(h, k2, k2);
                saxpy(0.5, k2, v, tmp);
                dV(t + h / 2, tmp, k3);
                sax(h, k3, k3);
                saxpy(1, k3, v, tmp);
                // Hack: Fake limit from right
                dV(s == steps - 1 ? t + h + 1e-14 : t + h, tmp, k4);
                sax(h, k4, k4);
                saxpy(1 / 6.0, k4, v, tmp);
                saxpy(2 / 6.0, k3, tmp, tmp);
                saxpy(2 / 6.0, k2, tmp, tmp);
                saxpy(1 / 6.0, k1, tmp, v);
            }
            Array.Copy(v, result[y - 1 - b], v.Length);
        }
        return result;
    }

    // sax = scalar a times x array, imperative version
    static void sax(double a, double[] x, double[] res)
    {
        if (x.Length != res.Length)
            throw new Exception("sax: lengths of x and res differ");
        for (int i = 0; i < x.Length; i++)
            res[i] = a * x[i];
    }

    // saxpy = scalar a times x array plus y array, imperative version
    static void saxpy(double a, double[] x, double[] y, double[] res)
    {
        if (x.Length != y.Length)
            throw new Exception("saxpy: lengths of x and y differ");
        if (x.Length != res.Length)
            throw new Exception("saxpy: lengths of x and res differ");
        for (int i = 0; i < x.Length; i++)
            res[i] = a * x[i] + y[i];
    }

    static void Print(double[][] result)
    {
        for (int y = 0; y < result.Length; y++)
        {
            Console.Write("{0,3}:", y);
            for (int i = 0; i < result[y].Length; i++)
                Console.Write("  {0,20:F16}", result[y][i]);
            //	Console.Write("  {0,20:F16} {1}", result[y][i], toBits(result[y][i]));
            Console.WriteLine();
        }
        Console.WriteLine();
    }

    // Need unsafe code to convert between 32-bit float and bits

    // unsafe public static String toBits(double d) {
    //   long* ptr = (long*) &d;
    //   StringBuilder sb = insertBlanks(toBits(*ptr, 64), 12, 1);
    //   return sb.ToString();
    // }

    public static String toBits(long n, int size)
    {
        char[] cs = new char[size];
        for (int i = 1; i <= size; i++)
        {
            cs[size - i] = (n & 1) != 0 ? '1' : '0';
            n >>= 1;
        }
        return new String(cs);
    }

    private static StringBuilder insertBlanks(String s, params int[] pos)
    {
        StringBuilder sb = new StringBuilder(s);
        for (int i = 0; i < pos.Length; i++)
            sb.Insert(pos[i], ' ');
        return sb;
    }


    // The two-state Actulus calculation kernel examples; 
    // really one-state because V1(t) = 0 for all t.

    class PureEndowment
    {
        static double b_0(double t)
        {
            return 0.0;
        }

        static double mu_01(double t)
        {
            return GM(t);
        }

        static double bj_00(double t)
        {
            // This works only because t is known to be an integer
            return t == pensiontime ? bpension : 0.0;
        }

        static double bj_01(double t)
        {
            return 0.0;
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n PureEndowment");
            // Print(new double[][] { new double[] { 0.14379469738 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
                 (double t, double[] res) => { res[0] = bj_00(t); },
                 40, 0, steps, new double[] { 0 });
        }
    }

    class DeferredTemporaryLifeAnnuity
    {
        static int m = 35, n = 10;

        static double b_0(double t)
        {
            return bpension * indicator(t > m) * indicator(t < m + n);
        }

        static double mu_01(double t)
        {
            return GM(t);
        }

        static double bj_00(double t)
        {
            return 0.0;
        }

        static double bj_01(double t)
        {
            return 0.0;
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n DeferredTemporaryLifeAnnuity");
            // Print(new double[][] { new double[] { 1.0265607675 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
                 (double t, double[] res) => { res[0] = bj_00(t); },
                 50, 0, steps, new double[] { 0 });
        }
    }

    class TemporaryLifeAnnuityPremium
    {
        static int n = 35;
        static double bpremium = 1;

        static double b_0(double t)
        {
            return -bpremium * indicator(t >= 0) * indicator(t < n);
        }

        static double mu_01(double t)
        {
            return GM(t);
        }

        static double bj_00(double t)
        {
            return 0.0;
        }

        static double bj_01(double t)
        {
            return 0.0;
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n TemporaryLifeAnnuityPremium");
            // Print(new double[][] { new double[] { -15.971767666 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
                 (double t, double[] res) => { res[0] = bj_00(t); },
                 50, 0, steps, new double[] { 0 });
        }
    }

    class TermInsurance
    {
        static int n = 35;
        static double bdeath = 1;

        static double b_0(double t)
        {
            return 0.0;
        }

        static double mu_01(double t)
        {
            return GM(t);
        }

        static double bj_00(double t)
        {
            return 0.0;
        }

        static double bj_01(double t)
        {
            return bdeath * indicator(t > 0) * indicator(t < n);
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n TermInsurance");
            // Print(new double[][] { new double[] { 0.057616919318 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            { res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t)); },
                 (double t, double[] res) => { res[0] = bj_00(t); },
                 50, 0, steps, new double[] { 0 });
        }
    }

    // The three-state Actulus calculation kernel examples; 
    // really two-state because V2(t) = 0 for all t.

    class DisabilityAnnuity
    {
        static int n = 35;
        static double bdisabled = 1;

        static double b_0(double t)
        {
            return 0.0;
        }

        static double b_1(double t)
        {
            return bdisabled * indicator(t > 0) * indicator(t < n);
        }

        const double ln10 = 2.3025850929940459;

        static double GM01(double t)
        {
            return 0.0006 + Math.Exp(ln10 * (4.71609 - 10 + 0.06 * (age + t)));
        }

        static double GM02(double t)
        {
            return GM(t);
        }

        static double GM12(double t)
        {
            return GM(t);
        }

        static double mu_01(double t)
        {
            return GM01(t);
        }

        static double mu_02(double t)
        {
            return GM02(t);
        }

        static double mu_12(double t)
        {
            return GM12(t);
        }

        static double bj_00(double t)
        {
            return 0.0;
        }

        static double bj_01(double t)
        {
            return 0.0;
        }

        static double bj_02(double t)
        {
            return 0.0;
        }

        static double bj_11(double t)
        {
            return 0.0;
        }

        static double bj_12(double t)
        {
            return 0.0;
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n DisabilityAnnuity");
            // Print(new double[][] { new double[] { 0.55552610797, 15.971767666 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            {
                res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t))
                         - mu_02(t) * (0 - V[0] + bj_02(t));
                res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t));
            },
                 (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
                 50, 0, steps, new double[] { 0, 0 });
        }
    }

    class DisabilityTermInsurance
    {
        static int n = 35;
        static double bdisabled = 1;

        static double b_0(double t)
        {
            return 0.0;
        }

        static double b_1(double t)
        {
            return 0.0;
        }

        const double ln10 = 2.3025850929940459;

        static double GM01(double t)
        {
            return 0.0006 + Math.Exp(ln10 * (4.71609 - 10 + 0.06 * (age + t)));
        }

        static double GM02(double t)
        {
            return GM(t);
        }

        static double GM12(double t)
        {
            return GM(t);
        }

        static double mu_01(double t)
        {
            return GM01(t);
        }

        static double mu_02(double t)
        {
            return GM02(t);
        }

        static double mu_12(double t)
        {
            return GM12(t);
        }

        static double bj_00(double t)
        {
            return 0.0;
        }

        static double bj_01(double t)
        {
            return bdisabled * indicator(t > 0) * indicator(t < n);
        }

        static double bj_02(double t)
        {
            return 0.0;
        }

        static double bj_11(double t)
        {
            return 0.0;
        }

        static double bj_12(double t)
        {
            return 0.0;
        }

        public static double[][] Compute()
        {
            // Console.WriteLine("\n DisabilityTermInsurance");
            // Print(new double[][] { new double[] { 0.071418699003, 0.000000000 } });
            // Console.WriteLine();
            return RK4_n((double t, double[] V, double[] res) =>
            {
                res[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t))
                         - mu_02(t) * (0 - V[0] + bj_02(t));
                res[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t));
            },
                 (double t, double[] res) => { res[0] = bj_00(t); res[1] = bj_11(t); },
                 50, 0, steps, new double[] { 0, 0 });
        }
    }

    public static void ComputeAll()
    {
        // Compute and print reserves
        Print(PureEndowment.Compute());
        Print(DeferredTemporaryLifeAnnuity.Compute());
        Print(TemporaryLifeAnnuityPremium.Compute());
        Print(TermInsurance.Compute());
        Print(DisabilityAnnuity.Compute());
        Print(DisabilityTermInsurance.Compute());
    }

    public static void ComputeAllNoPrint()
    {
        // Compute reserves
        PureEndowment.Compute();
        DeferredTemporaryLifeAnnuity.Compute();
        TemporaryLifeAnnuityPremium.Compute();
        TermInsurance.Compute();
        DisabilityAnnuity.Compute();
        DisabilityTermInsurance.Compute();
    }

    public static void Main(String[] args)
    {
        // double[] result = RK4_1(40, 0, 256, 0);
        // for (int year=0; year<=40; year++)
        //   Console.WriteLine("{0,3}: {1,20:F16}", year, result[year]);
        // Print(RK4_2(40, 0, 256, 0, 0));
        // Print(RK4_n(40, 0, 256, new double[] { 0, 0 }));
        // for (double year=0; year<=40; year+=0.25)
        //   Console.WriteLine("{0,7:F2}: {1,10:F7}", year, rFsa(year));
        CalculationSpecifications.ComputeAll();
        Timer t = new Timer();
        const int count = 1000;
        for (int i = 0; i < count; i++)
            CalculationSpecifications.ComputeAllNoPrint();
        double time = t.Check();
        Console.WriteLine("{0,8:F3} ms", 1000.0 * time / count);
    }
}

public class Timer
{
    private Stopwatch stopwatch;

    public Timer()
    {
        stopwatch = new Stopwatch();
        stopwatch.Reset();
        stopwatch.Start();
    }

    public double Check()
    {
        return stopwatch.ElapsedMilliseconds / 1000.0;
    }
}
