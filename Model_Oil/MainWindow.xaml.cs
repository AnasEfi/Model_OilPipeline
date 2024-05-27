using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using static MaterialDesignThemes.Wpf.Theme;

namespace Model_Oil
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public static MainWindow Instance { get; private set; }
        public MainWindow()
        {
            InitializeComponent();
            Instance = this;
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {

            double nu_20 = double.Parse(nu_20_text.Text.Trim());
            double nu_50 = double.Parse(nu_50_text.Text.Trim());
            double po_20  = double.Parse(po_20_text.Text.Trim());
            double D_out = double.Parse(D_out_text.Text.Trim()); 
            double thick = double.Parse(thick_text.Text.Trim()); 
            double delta = double.Parse(delta_text.Text.Trim()); 
            double L = double.Parse(L_text.Text.Trim());
            double K_mn = double.Parse(K_mn_text.Text.Trim());
            double eps = double.Parse(eps_text.Text.Trim());
            double T_soil = double.Parse(T_soil_text.Text.Trim());
            double P_end = double.Parse(P_end_text.Text.Trim());
            double T_start = double.Parse(T_start_text.Text.Trim());
            double Q =  double.Parse(Q_text.Text.Trim());
            
            ModelEquation modelNew = new ModelEquation(nu_20, nu_50, po_20, D_out,
             thick, delta, L, K_mn, eps, T_soil , P_end, T_start, Q);

            modelNew.BeginCalculations();


        }
    }
}