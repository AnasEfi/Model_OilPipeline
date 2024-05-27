using System;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Threading.Tasks;

public class ModelEquation
{
    public const double T_20 = 20.0 + 273;
    public const double T_50 = 50.0 + 273;
    public const double tolerance_T = 0.000001;
    public const double tolerance_p = 0.00001;
    public const double g = 9.81;
    public const double delta_Z = 0.0;

    double nu_20;
    double nu_50;
    double po_20;
    double D_out;
    double thick;
    double delta;
    double L;
    double K_mn;
    double eps;
    double T_soil;
    double P_start;
    double P_end;
    double T_start;
    double Q;

    public ModelEquation(double nu_20, double nu_50, double po_20, double D_out,
        double thick, double delta, double L, double K_mn, double eps, double T_soil, double P_start,
        double P_end, double T_start, double Q)
    {
        this.nu_20 = nu_20 * Math.Pow(10, -6);
        this.nu_50 = nu_50 * Math.Pow(10, -6);
        this.po_20 = po_20;
        this.D_out = D_out * 0.001;
        this.thick = thick * 0.001;
        this.delta = delta * 0.001;
        this.L = L * 100;
        this.K_mn = K_mn;
        this.eps = eps;
        this.T_soil = T_soil + 273;
        this.P_start = P_start;
        this.P_end = P_end * Math.Pow(10, 6);
        this.T_start = T_start + 273;
        this.Q = Q / 3600;
    }

    public void BeginCalculations()
    {
        double G = 0.0;
        double v_avg = 0.0;
        double P_start = 0.0;
        double H_in = 0.0;
        double H_out = 0.0;
        double hydroclone = 0.0;

        double Re_imper = 8.0 * Math.Pow(10, -5) * Math.Pow(D_out * 1000, 3) - 0.304 * Math.Pow(D_out * 1000, 2) + 412.95 * D_out * 1000 - 66028;
        double A = -3.0 * Math.Pow(10, -12) * Math.Pow(D_out * 1000, 3) + Math.Pow(10, -8) * Math.Pow(D_out * 1000, 2) - Math.Pow(10, -5) * D_out * 1000 + 0.0166;
        double k_t = KtCoeffCalc();
        double D_int = D_out - 2.0 * thick;

        double T_end = 0;
        double viscosity = 0.0;
        double Re = 0.0;
        double L_n = 0.0;
        double lambda = 0.0;
        double delta_P = 0.0;
        double delta_P_new = 0.0;
        double aL = 0.0;
        double c_v = 0.0;
        double po_avg = 0.0;
        double viscosity_start = 0.0;
        double T_end_new;
        int iter = 0;


        double T_avg = T_start;

        for (; iter < 100000; iter++)
        {
            po_avg = poFunc(T_avg);
            c_v = c_vCalc(T_avg);
            viscosity = viscosityFunc(T_avg, k_t);
            Re = ReCalc(viscosity, D_int);
            lambda = lambdaCalc(Re_imper, Re, A, iter);
            L_n = L + eps * D_int / lambda;   // Эквивалентная (расчетная) длина участка трубопровода
            if (iter == 0)
            {
                delta_P = delta_PCalc(lambda, D_int, L_n, po_avg); // Расчет разницы давлений на входе и выходе
                aL = CalcShukhovCriterion_aL(po_avg, c_v);
                T_avg = T_avgCalc(aL, c_v, delta_P, po_avg);
                T_end = T_finalCalc(aL, po_avg, c_v, delta_P_new);
            }
            else
            {
                delta_P_new = delta_PCalc(lambda, D_int, L_n, po_avg); // Расчет разницы давлений на входе и выходе
                T_end_new = T_finalCalc(aL, po_avg, c_v, delta_P_new);

                if (ContinueCondition(delta_P, delta_P_new))
                {
                    delta_P = delta_P_new;
                    T_end = T_end_new;
                }
                else
                {
                    G = GCalc(po_avg);
                    v_avg = OilSpeedCalc(G, po_avg, D_int);
                    P_start = P_end + delta_P_new;
                    H_in = HCalc(po_avg, P_start);
                    H_out = HCalc(po_avg, P_end);
                    hydroclone = HydrocloneCalc(lambda, D_int, v_avg);
                    T_end = T_finalCalc(aL, po_avg, c_v, delta_P_new);
                    T_avg = T_avgCalc(aL, c_v, delta_P, po_avg);
                    viscosity_start = viscosityFunc(T_start, k_t);
                    Console.WriteLine($"The end of iterations, iter : {iter}");
                    break;
                }
            }
        }
        Console.WriteLine($"Found params: \n" +
            $"T_start = {T_start}\nT_out = {T_end},\nT_avg = {T_avg},\npo_avg = {po_avg},\nc_avg = {c_v}\n" +
            $"v_avg = {v_avg},\nRe = {Re},\nlambda = {lambda},\niter = {iter},\n" +
            $"G = {G},\nP_start = {P_start},\nP_out = {P_end},\nviscosity_start = {viscosity_start} \n" +
            $"viscosity_end = {viscosity},\nH_in = {H_in},\nH_out = {H_out},\n" +
            $"Hydroclone = {hydroclone},\nT_avg = {T_avg},\nT_final = {T_end}");
    }
    public double KtCoeffCalc()
    {
        return Math.Log(nu_20 / nu_50) / (T_50 - T_20);
    }
    public double T_finalCalc(double aL, double po_avg, double c_v, double delta_P)
    {
        double T_final = T_soil + (T_start - T_soil) * Math.Exp(-aL) + (1 * delta_P * Math.Pow(10, -6)) / (po_avg * c_v * aL)
            * (1 - Math.Exp(-aL)) - g * delta_Z / (c_v * aL) * (1 - Math.Exp(-aL));
        return T_final;
    }
    public double poFunc(double T)
    {
        return po_20 + ((1.831 - 0.0013233 * po_20) * (293.15 - T));
    }
    public double viscosityFunc(double T, double k_t)
    {
        //Формула Филонова-Рейнольдса
        return nu_20 * Math.Exp(-k_t * (T - 293.15));
    }
    public double lambdaCalc(double Re_imper, double Re, double A, double iter)
    {
        Console.WriteLine("Now iter is : " + iter);
        if (Re < 2000)
            return 64 / Re;
        else if (2000 < Re && Re < 2800)
            return (0.16 * Re - 13) * Math.Pow(10.0, -4);
        else if (2800 < Re && Re < Re_imper)
            return (0.3164 / Math.Pow(Re, 1.0 / 4.0));
        else if (Re > Re_imper)
            return A + 1.7 / Math.Sqrt(Re);
        else
        {
            Console.WriteLine("The val Re is: " + Re);
            throw new Exception("Error param Re");
        }
    }

    // Удельная изохорная теплоемкость
    public double c_vCalc(double T)
    {
        Console.WriteLine("C_v is: " + 31.56 / (Math.Sqrt(po_20)) * (762 + 3.93 * T));
        return 31.56 / (Math.Sqrt(po_20)) * (762 + 3.93 * T);
    }
    public double delta_PCalc(double lambda, double D_int, double L_n, double po_avg)
    {
        return lambda * 8.0 * Math.Pow(Q, 2) / (Math.Pow(Math.PI, 2) * Math.Pow(D_int, 5)) * po_avg * L_n + po_avg * delta_Z * g; // delta Z = 0
    }
    public double ReCalc(double nu_avg, double D_int)
    {
        return 4.0 * Q / (3.14 * D_int * nu_avg);
    }
    public double CalcShukhovCriterion_aL(double po_avg, double c_v)
    {
        return ((Math.PI * K_mn * D_out) / (po_avg * Q * c_v)) * L;
    }
    public double T_avgCalc(double aL, double c_v, double delta_P, double po_avg)
    {
        return T_soil + (T_start - T_soil) * ((1.0 - Math.Exp(-aL)) / aL)
            + g / (c_v * aL) * (delta_P * Math.Pow(10, -6) / (po_avg * g) - delta_Z) * (1.0 - ((1.0 - Math.Exp(-aL)) / aL));
    }

    public bool ContinueCondition(double delta_P, double delta_P_new)
    {
        if (((Math.Abs(delta_P_new - delta_P) / delta_P_new)) > tolerance_p)
            return true;
        else return false;
    }

    public double GCalc(double po_avg)
    {
        return po_avg * Q;
    }
    public double OilSpeedCalc(double G, double po_avg, double D_int)
    {
        return 4.0 * G / (po_avg * Math.PI * Math.Pow(D_int, 2));
    }
    public double HCalc(double po_avg, double p_in_out)
    {
        return p_in_out / (po_avg * g);
    }
    public double HydrocloneCalc(double lambda, double D_int, double v_avg)
    {
        return lambda * Math.Pow(v_avg, 2) / (D_int * 2 * g);
    }

}
