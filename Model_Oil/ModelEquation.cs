using System;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Threading.Tasks;

namespace Model_Oil
{
    public class ModelEquation
    {
        public const double T_20 = 20.0 + 273;
        public const double T_50 = 50.0 + 273;
        public const double tolerance_p = 0.000001;
        public const double tolerance_t = 0.000001;

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
            double thick, double delta, double L, double K_mn, double eps, double T_soil,
            double P_end, double T_start, double Q)
        {
            this.nu_20 = nu_20 * Math.Pow(10, -6);
            this.nu_50 = nu_50 * Math.Pow(10, -6);
            this.po_20 = po_20;
            this.D_out = D_out * 0.001;
            this.thick = thick * 0.001;
            this.delta = delta * 0.001;
            this.L = L * 1000;
            this.K_mn = K_mn;         //коэффициент теплопередачи
            this.eps = eps;
            this.T_soil = T_soil + 273;
            //this.P_start = P_start;
            this.P_end = P_end * Math.Pow(10, 6);
            this.T_start = T_start + 273;
            this.Q = Q / 3600;
        }

        public void BeginCalculations()
        {

            double G = 0.0;
            double speed_avg = 0.0;
            double P_start = 0.0;
            double H_in = 0.0;
            double H_out = 0.0;
            double hydroclone = 0.0;
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
            int iter = 0;
            double T_avg_new = 0;

            double Re_imper = 8.0 * Math.Pow(10, -5) * Math.Pow(D_out * 1000, 3) - 0.304 * Math.Pow(D_out * 1000, 2) + 412.95 * D_out * 1000 - 66028;
            double A = -3.0 * Math.Pow(10, -12) * Math.Pow(D_out * 1000, 3) + Math.Pow(10, -8) * Math.Pow(D_out * 1000, 2) - Math.Pow(10, -5) * D_out * 1000 + 0.0166;
            double k_t = KtCoeffCalc(); 
            double D_int = D_out - 2.0 * thick;
            double T_avg = T_start;

            for (; iter < 1000; iter++)
            {
                po_avg = poFunc(T_avg);
                c_v = c_vCalc(T_avg);
                viscosity = viscosityFunc(T_avg, k_t); // кинематическая вязкость 
                Re = ReCalc(viscosity, D_int); 
                lambda = lambdaCalc(Re_imper, Re, A); // коэффициент гидравлического сопростивления
                L_n = L + eps * D_int / lambda;   // Эквивалентная (расчетная) длина участка трубопровода
                if (iter == 0)
                {
                    delta_P = delta_PCalc(lambda, D_int, L_n, po_avg); // Расчет разницы давлений на входе и выходе
                    aL = CalcShukhovCriterion_aL(po_avg, c_v);
                    T_avg = T_avgCalc(aL, c_v, delta_P, po_avg);
                }
                else
                {
                    delta_P_new = delta_PCalc(lambda, D_int, L_n, po_avg); // Расчет разницы давлений на входе и выходе
                    T_avg_new = T_avgCalc(aL, c_v, delta_P_new, po_avg);

                    if (ContinueCondition(delta_P, delta_P_new, T_avg, T_avg_new))
                    {
                        delta_P = delta_P_new;
                        T_avg = T_avg_new;
                    }
                    else
                    {
                        G = GCalc(po_avg);
                        speed_avg = OilSpeedCalc(G, po_avg, D_int);
                        P_start = P_end + delta_P_new;
                        H_in = HCalc(poFunc(T_start), P_start);
                        H_out = HCalc(poFunc(T_end), P_end);
                        hydroclone = HydrocloneCalc(lambda, D_int, speed_avg);
                        T_end = T_finalCalc(aL, po_avg, c_v, delta_P_new);
                        T_avg = T_avgCalc(aL, c_v, delta_P, po_avg);
                        break;
                    }
                }
            }
            
            if (iter < 1000)
            {
                MainWindow.Instance.T_avg_text.Text = Math.Round(T_avg, 2).ToString();
                MainWindow.Instance.po_avg_text.Text = Math.Round(po_avg, 2).ToString();
                MainWindow.Instance.speed_avg_text.Text = Math.Round(speed_avg, 2).ToString();

                var visc = viscosityFunc(T_avg, k_t);
                MainWindow.Instance.nu_avg_text.Text = Math.Round(visc, 10).ToString();
                MainWindow.Instance.c_v_text.Text = Math.Round(c_vCalc(T_avg), 2).ToString();
                MainWindow.Instance.Re_text.Text = Math.Round(Re, 2).ToString();

                MainWindow.Instance.P_start_text_2.Text = Math.Round(P_start * Math.Pow(10, -6), 2).ToString();
                MainWindow.Instance.P_end_text_2.Text = Math.Round(P_end * Math.Pow(10,-6), 2).ToString();
                MainWindow.Instance.T_start_text_2.Text = Math.Round(T_start - 273, 2).ToString();
                MainWindow.Instance.T_end_text.Text = Math.Round(T_end-273, 2).ToString();
                MainWindow.Instance.po_start_text.Text = Math.Round(poFunc(T_start), 2).ToString();

                MainWindow.Instance.po_end_text.Text = Math.Round(poFunc(T_end), 2).ToString();
                MainWindow.Instance.nu_start_text.Text = Math.Round(viscosityFunc(T_start, k_t) * Math.Pow(10,4), 6).ToString();
                MainWindow.Instance.nu_end_text.Text = Math.Round(viscosityFunc(T_end, k_t) * Math.Pow(10, 4), 6).ToString();
                MainWindow.Instance.speed_start_text.Text = Math.Round(OilSpeedCalc(G, poFunc(T_start), D_int), 6).ToString();
                MainWindow.Instance.speed_end_text.Text = Math.Round(OilSpeedCalc(G, poFunc(T_end), D_int), 6).ToString();
                MainWindow.Instance.H_start_text.Text = Math.Round(HCalc(poFunc(T_start), P_start), 2).ToString();

                MainWindow.Instance.H_end_text.Text = Math.Round(HCalc(poFunc(T_end), P_end), 2).ToString();
                MainWindow.Instance.Q_text_2.Text = MainWindow.Instance.Q_text.Text;
                MainWindow.Instance.G_text.Text = Math.Round(G * 2 / 1000 * 3600, 2).ToString();
                MainWindow.Instance.DiffH_text.Text = (H_in-H_out).ToString();
                MainWindow.Instance.Hydroclone_text.Text = Math.Round(hydroclone, 10).ToString();
                MainWindow.Instance.iter_text.Text = iter.ToString();
                
            }
            else
            {
                MainWindow.Instance.iter_text.Text = "Превышен лимит итераций".ToString();
            }
        }
        public double KtCoeffCalc()
        {
            return Math.Log(nu_20 / nu_50) / ((T_50-273) - (T_20-273));
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
        public double lambdaCalc(double Re_imper, double Re, double A)
        {
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
                throw new Exception("Error param Re");
            }
        }

        // Удельная изохорная теплоемкость
        public double c_vCalc(double T)
        {
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

        public bool ContinueCondition(double delta_P, double delta_P_new, double T_avg, double T_avg_new)
        {
            if (((Math.Abs(delta_P - delta_P_new) / delta_P) > tolerance_p) || ((Math.Abs(T_avg - T_avg_new) / T_avg) > tolerance_t))
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
            return lambda * (1 / D_int) * (Math.Pow(v_avg, 2) / (2 * g)) * 1000; 
        }
    }
}
