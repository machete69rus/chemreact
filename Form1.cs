using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Globalization;
using System.Text.RegularExpressions;

namespace WindowsFormsApp1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        public class CompoundInfo
        { 
            public string Formula { get; set; }
            public int Coefficient { get; set; }

            public double MolarMass { get; set; }
        }

        private List<CompoundInfo> lastBalancedCompounds = new List<CompoundInfo>();

        // Словарь с атомными массами химических элементов, взято из
        // https://ru.wikipedia.org/wiki/Список_химических_элементов
        private readonly Dictionary<string, double> PeriodicTable = new Dictionary<string, double>() 
        {
            // Период 1
            { "H", 1.00794 },{ "He", 4.002602 },
            
            // Период 2
            { "Li", 6.941 }, { "Be", 9.012182 }, { "B", 10.811 }, { "C", 12.0107 }, 
            { "N", 14.0067 }, { "O", 15.9994 }, { "F", 18.9984032 }, { "Ne", 20.1797 },

            // Период 3
            { "Na", 22.98976928 }, { "Mg", 24.3050 }, { "Al", 26.9815386 }, { "Si", 28.0855 },
            { "P", 30.973762 }, { "S", 32.065 }, { "Cl", 35.453 }, { "Ar", 39.948 },

            // Период 4
            { "K", 39.0983 }, { "Ca", 40.078 }, { "Sc", 44.955912 }, { "Ti", 47.867 },
            { "V", 50.9415 }, { "Cr", 51.9961 }, { "Mn", 54.938045 }, { "Fe", 55.845 },
            { "Co", 58.933195 }, { "Ni", 58.6934 }, { "Cu", 63.546 }, { "Zn", 65.409 },
            { "Ga", 69.723 }, { "Ge", 72.64 }, { "As", 74.92160 }, { "Se", 78.96 },
            { "Br", 79.904 }, { "Kr", 83.798 },

            // Период 5
            { "Rb", 85.4678 }, { "Sr", 87.62 }, { "Y", 88.90585 }, { "Zr", 91.224 },
            { "Nb", 92.90638 }, { "Mo", 95.94 }, { "Tc", 98.9063 }, { "Ru", 101.07 },
            { "Rh", 102.90550 }, { "Pd", 106.42 }, { "Ag", 107.8682 }, { "Cd", 112.411 },
            { "In", 114.818 }, { "Sn", 118.710 }, { "Sb", 121.760 }, { "Te", 127.60 },
            { "I", 126.90447 }, { "Xe", 131.293 },

            // Период 6

            { "Cs", 132.9054519 }, { "Ba", 137.327 },

                 // Лантаноиды
                { "La", 138.90547 }, { "Ce", 140.116 }, { "Pr", 140.90765 }, { "Nd", 144.242 },
                { "Pm", 146.9151 }, { "Sm", 150.36 }, { "Eu", 151.964 }, { "Gd", 157.25 },
                { "Tb", 158.92535 }, { "Dy", 162.500 }, { "Ho", 164.93032 }, { "Er", 167.259 },
                { "Tm", 168.93421 }, { "Yb", 173.04 }, { "Lu", 174.967 },

            { "Hf", 178.49 }, { "Ta", 180.9479 }, { "W", 183.84 }, { "Re", 186.207 },
            { "Os", 190.23 }, { "Ir", 192.217 }, { "Pt", 195.084 }, { "Au", 196.966569 },
            { "Hg", 200.59 }, { "Tl", 204.3833 }, { "Pb", 207.2 }, { "Bi", 208.98040 },
            { "Po", 208.9824 }, { "At", 209.9871 }, { "Rn", 222.0176 }, 
            
            // Период 7
            { "Fr", 223.0197 }, { "Ra", 226.0254 },

                //Актиноиды
                { "Ac", 227.0278 }, { "Th", 232.03806 }, { "Pa", 231.03588 }, { "U", 238.02891 },
                { "Np", 237.0482 }, { "Pu", 244.0642 }, { "Am", 243.0614 }, { "Cm", 247.0703 },
                { "Bk", 247.0703 }, { "Cf", 251.0796 }, { "Es", 252.0829 }, { "Fm", 257.0951 },
                { "Md", 258.0986 }, { "No", 259.1009 },

            { "Lr", 266.00 }, { "Rf", 267.00 }, { "Db", 268.00 }, { "Sg", 269.00 },
            { "Bh", 270.00 }, { "Hs", 277.00 }, { "Mt", 278.00 }, { "Ds", 281.00 },
            { "Rg", 282.00 }, { "Cn", 285.00 }, { "Nh", 286.00 }, { "Fl", 289.00 },
            { "Mc", 290.00 }, { "Lv", 293.00 }, { "Ts", 294.00 }, { "Og", 294.00 },

        };

        Dictionary<string, double> reagentPurity = new Dictionary<string, double>();

        private List<string> allCompounds;
        public List<int> balanceCoefficients;
        private string[] reactants, products;

        private Dictionary<string, (int Coefficient, double MolarMass)> compoundData = new Dictionary<string, (int, double)>();

        private void balanceButton_Click(object sender, EventArgs e)
        {
            try
            {
                var (result, log) = BalanceEquation(inputBox.Text);
                outputBox.Text = result;
                ShowLog(log);
               
            }
            catch (Exception ex)
            {
                MessageBox.Show("Ошибка: " + ex.Message);
            }
        }

        /*
        private string BalanceEquation(string equation)
        {
            equation = equation.Replace("→", "=>").Replace(" ", "");
            var sides = equation.Split(new[] { "=>", "=" }, StringSplitOptions.RemoveEmptyEntries);
            if (sides.Length != 2) throw new Exception("Уравнение должно содержать '=>'");

            var reactants = sides[0].Split('+');
            var products = sides[1].Split('+');
            var compounds = reactants.Concat(products).ToList();

            var elements = new HashSet<string>();
            foreach (var compound in compounds)
                foreach (var (el, _) in ParseCompound(compound))
                    elements.Add(el);

            var elementList = elements.ToList();
            int rows = elementList.Count;
            int cols = compounds.Count;

            double[,] matrix = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                string element = elementList[i];

                for (int j = 0; j < reactants.Length; j++)
                    matrix[i, j] = GetElementCount(reactants[j], element);

                for (int j = 0; j < products.Length; j++)
                    matrix[i, j + reactants.Length] = -GetElementCount(products[j], element);
            }

            var solution = SolveHomogeneousSystem(matrix);
            var intSolution = ScaleToIntegers(solution);

            if (intSolution.Count != compounds.Count)
                throw new Exception("Ошибка при решении уравнения");

            string left = string.Join(" + ", reactants.Select((r, i) => (intSolution[i] == 1 ? "" : intSolution[i] + "") + r));
            string right = string.Join(" + ", products.Select((p, i) => (intSolution[i + reactants.Length] == 1 ? "" : intSolution[i + reactants.Length] + "") + p));

           

            return left + " => " + right;
        }

        */

        private (string, List<(string message, bool isBalanced)>) BalanceEquation(string equation)
        {
            equation = equation.Replace("→", "=>").Replace(" ", "");
            var sides = equation.Split(new[] { "=>", "=" }, StringSplitOptions.RemoveEmptyEntries);
            if (sides.Length != 2) throw new Exception("Уравнение должно содержать '=>'");

            var reactants = sides[0].Split('+');
            var products = sides[1].Split('+');
            var compounds = reactants.Concat(products).ToList();

            var elements = new HashSet<string>();
            foreach (var compound in compounds)
                foreach (var (el, _) in ParseCompound(compound))
                    elements.Add(el);

            var elementList = elements.ToList();
            int rows = elementList.Count;
            int cols = compounds.Count;

            double[,] matrix = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                string element = elementList[i];

                for (int j = 0; j < reactants.Length; j++)
                    matrix[i, j] = GetElementCount(reactants[j], element);

                for (int j = 0; j < products.Length; j++)
                    matrix[i, j + reactants.Length] = -GetElementCount(products[j], element);
            }

            var solution = SolveHomogeneousSystem(matrix);
            var intSolution = ScaleToIntegers(solution);

            if (intSolution.Count != compounds.Count)
                throw new Exception("Ошибка при решении уравнения");

            string left = string.Join(" + ", reactants.Select((r, i) => (intSolution[i] == 1 ? "" : intSolution[i] + "") + r));
            string right = string.Join(" + ", products.Select((p, i) => (intSolution[i + reactants.Length] == 1 ? "" : intSolution[i + reactants.Length] + "") + p));

            // Формируем лог по каждому элементу
            var log = new List<(string, bool)>();
            foreach (var element in elementList)
            {
                double leftCount = 0, rightCount = 0;

                for (int i = 0; i < reactants.Length; i++)
                    leftCount += GetElementCount(reactants[i], element) * intSolution[i];

                for (int i = 0; i < products.Length; i++)
                    rightCount += GetElementCount(products[i], element) * intSolution[i + reactants.Length];

                bool balanced = Math.Abs(leftCount - rightCount) < 1e-6;
                log.Add((
                    $"{element}: слева {leftCount}, справа {rightCount} — " + (balanced ? "сбалансировано" : "НЕСБАЛАНСИРОВАНО"),
                    balanced));
            }
            this.reactants = reactants;
            this.products = products;
            this.allCompounds = reactants.Concat(products).ToList();
            this.balanceCoefficients = intSolution;

            dataGridView1.Rows.Clear();
            foreach (var compound in reactants)
            {
                dataGridView1.Rows.Add(compound, "100");
            }

            /*
            lastBalancedCompounds = compounds.Select((compound, i) => new CompoundInfo
                {
                    Formula=compound, Coefficient = intSolution[i], MolarMass = CalculateMolarMass(compound)
                }
                    ).ToList();
            */

            // Рабочий кусок
            /*
            productscomboBox.Items.Clear();
           
            foreach (var compound in lastBalancedCompounds)
                productscomboBox.Items.Add(compound.Formula);

            UpdateProductComboBox();
            */

            lastBalancedCompounds.Clear();
            productscomboBox.Items.Clear();

            string Left = "";
            string Right = "";

            for (int i = 0; i < reactants.Length; i++)
            {
                int coeff = intSolution[i];
                string formula = reactants[i];

                Left += (coeff == 1 ? "" : coeff + "") + formula;
                if (i < reactants.Length - 1) Left += " + ";

                lastBalancedCompounds.Add(new CompoundInfo
                {
                    Formula = formula,
                    Coefficient = coeff,
                    MolarMass = CalculateMolarMass(formula)
                });
                productscomboBox.Items.Add(formula);
            }

            for (int i = 0; i < products.Length; i++)
            {
                int coeff = intSolution[i + reactants.Length];
                string formula = products[i];

                Right += (coeff == 1 ? "" : coeff + "") + formula;
                if (i < products.Length - 1) Right += " + ";

                lastBalancedCompounds.Add(new CompoundInfo
                {
                    Formula = formula,
                    Coefficient = coeff,
                    MolarMass = CalculateMolarMass(formula)
                });
                productscomboBox.Items.Add(formula);
            }
            
            compoundData.Clear();

            for (int i = 0; i < compounds.Count; i++)
            {
                string compound = compounds[i];
                int coefficient = intSolution[i];
                double molarMass = CalculateMolarMass(compound);

                compoundData[compound] = (coefficient, molarMass);
            }

            return (left + " => " + right, log);
           
        }

        private double CalculateMolarMass(string formula)
        {
            var atoms = ParseCompound(formula); // ты уже используешь ParseCompound
            double mass = 0;

            foreach (var (element, count) in atoms)
            {
                if (!PeriodicTable.ContainsKey(element))
                    throw new Exception($"Неизвестный элемент: {element}");

                mass += PeriodicTable[element] * count;
            }

            return mass;
        }

        private void UpdateProductComboBox()
        {
            if (products == null || balanceCoefficients == null || products.Length == 0)
                return;

            productscomboBox.Items.Clear();

            for (int i = 0; i < products.Length; i++)
            {
                int globalIndex = i + reactants.Length;
                if (globalIndex < balanceCoefficients.Count)
                {
                    string itemText = $"{products[i]} (коэф. {balanceCoefficients[globalIndex]})";
                    productscomboBox.Items.Add(itemText);
                }
            }

            if (productscomboBox.Items.Count > 0)
                productscomboBox.SelectedIndex = 0;
        }

        private List<(string element, double count)> ParseCompound(string compound)
        {
            var list = new List<(string, double)>();
            var matches = Regex.Matches(compound, @"([A-Z][a-z]*)(\d*\.?\d*)");

            foreach (Match match in matches)
            {
                string el = match.Groups[1].Value;
                string countStr = match.Groups[2].Value;
                double count = 1;

                if (!string.IsNullOrEmpty(countStr))
                    count = double.Parse(countStr, CultureInfo.InvariantCulture);

                list.Add((el, count));
            }

            return list;
        }
       
        private double GetElementCount(string compound, string element)
        {
            return ParseCompound(compound).Where(x => x.element == element).Sum(x => x.count);
        }

        private double[] SolveHomogeneousSystem(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            double[,] augmented = new double[rows, cols + 1];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    augmented[i, j] = matrix[i, j];

            for (int i = 0; i < rows; i++)
                augmented[i, cols] = 0; // Right-hand side is 0 for homogeneous system

            int rank = GaussianElimination(augmented, rows, cols);

            double[] solution = new double[cols];
            for (int i = 0; i < cols; i++)
                solution[i] = 1;

            for (int i = rank - 1; i >= 0; i--)
            {
                int pivotCol = -1;
                for (int j = 0; j < cols; j++)
                {
                    if (Math.Abs(augmented[i, j]) > 1e-8)
                    {
                        pivotCol = j;
                        break;
                    }
                }

                if (pivotCol == -1) continue;

                double sum = 0;
                for (int j = pivotCol + 1; j < cols; j++)
                    sum += augmented[i, j] * solution[j];

                solution[pivotCol] = -sum / augmented[i, pivotCol];
            }

            return solution;
        }

        private int GaussianElimination(double[,] matrix, int rows, int cols)
        {
            int rank = 0;

            for (int col = 0; col < cols && rank < rows; col++)
            {
                int pivot = rank;
                for (int i = rank + 1; i < rows; i++)
                    if (Math.Abs(matrix[i, col]) > Math.Abs(matrix[pivot, col]))
                        pivot = i;

                if (Math.Abs(matrix[pivot, col]) < 1e-8) continue;

                SwapRows(matrix, rank, pivot, cols + 1);

                double div = matrix[rank, col];
                for (int j = 0; j <= cols; j++)
                    matrix[rank, j] /= div;

                for (int i = 0; i < rows; i++)
                {
                    if (i != rank && Math.Abs(matrix[i, col]) > 1e-8)
                    {
                        double factor = matrix[i, col];
                        for (int j = 0; j <= cols; j++)
                            matrix[i, j] -= factor * matrix[rank, j];
                    }
                }

                rank++;
            }

            return rank;
        }

        private void SwapRows(double[,] matrix, int row1, int row2, int totalCols)
        {
            for (int j = 0; j < totalCols; j++)
            {
                double temp = matrix[row1, j];
                matrix[row1, j] = matrix[row2, j];
                matrix[row2, j] = temp;
            }
        }

        private List<int> ScaleToIntegers(double[] values)
        {
            const int maxDenominator = 1000; // повысим точность

            var scaled = values.Select(v =>
            {
                // Умножаем на 1000, округляем, делим обратно
                double rounded = Math.Round(v * maxDenominator);
                return rounded;
            }).ToList();

            long gcd = scaled.Select(v => (long)Math.Abs(v)).Where(v => v > 0).Aggregate(GCD);

            if (gcd == 0) gcd = 1;

            return scaled.Select(v => (int)(v / gcd)).ToList();
        }

        /*
        private List<int> ScaleToIntegers(double[] values)
        {
            const int precision = 1000000;

            var scaled = values.Select(v => Math.Round(v * precision)).ToList();
            var gcd = scaled.Select(v => (long)Math.Abs(v)).Aggregate(GCD);

            if (gcd == 0) gcd = 1;

            return scaled.Select(v => (int)(v / gcd)).ToList();
        }
        */

        private long GCD(long a, long b)
        {
            return b == 0 ? a : GCD(b, a % b);
        }

        private void ShowLog(List<(string message, bool isBalanced)> log)
        {
            logBox.Clear();
            foreach (var (message, isBalanced) in log)
            {
                logBox.SelectionColor = isBalanced ? Color.DarkGreen : Color.Red;
                logBox.AppendText(message + Environment.NewLine);
            }
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
           
        }

        private double GetMolarMass(string formula)  // тут вторая таблица ПСХЭ
        {
            var periodicTable = new Dictionary<string, double>
    {
        { "H", 1.008 }, { "O", 16.00 }, { "Fe", 55.845 }, { "C", 12.01 },
        { "N", 14.007 }, { "S", 32.06 }, { "Na", 22.99 }, { "Cl", 35.45 },
        // Добавляй по мере необходимости
    };

            double molarMass = 0;
            foreach (var (el, count) in ParseCompound(formula))
            {
                if (!periodicTable.ContainsKey(el))
                    throw new Exception($"Неизвестный элемент: {el}");

                molarMass += periodicTable[el] * count;
            }

            return molarMass;
        }

        private string CalculateRequiredMasses(string selectedProduct, double desiredMass)
        {
            int index = products.ToList().FindIndex(p => selectedProduct.StartsWith(p));
            if (index == -1) throw new Exception("Продукт не найден");

            int productGlobalIndex = index + reactants.Length;
            double productCoeff = balanceCoefficients[productGlobalIndex];
            double productMolarMass = GetMolarMass(products[index]);

            double desiredMoles = desiredMass / productMolarMass;
            double moleFactor = desiredMoles / productCoeff;

            StringBuilder log = new StringBuilder();
            log.AppendLine($"Цель: получить {desiredMass} г {products[index]} ({desiredMoles:0.###} моль)");

            for (int i = 0; i < reactants.Length; i++)
            {
                string reactant = reactants[i];
                int coeff = balanceCoefficients[i];
                double moles = coeff * moleFactor;
                double mass = moles * GetMolarMass(reactant);
                log.AppendLine($"{reactant}: {mass:0.###} г (нужно {moles:0.###} моль)");
            }

            return log.ToString();
        }

        private void UpdateReagentPurities()
        {
            reagentPurity.Clear();
            foreach (var kvp in compoundData)
            {
                reagentPurity[kvp.Key] = 100.0; // По умолчанию — 100%
            }
        }

        private void checkBoxPurity_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Enabled = checkBoxPurity.Checked;
        }
        private void CalculateMassWithPurity()
        {
            if (productscomboBox.SelectedItem == null || string.IsNullOrWhiteSpace(massInputBox.Text))
            {
                MessageBox.Show("Выберите вещество и введите массу.");
                return;
            }

            string targetCompound = productscomboBox.SelectedItem.ToString();
            if (!compoundData.ContainsKey(targetCompound))
            {
                MessageBox.Show("Выбранное вещество не найдено в уравнении.");
                return;
            }

            if (!double.TryParse(massInputBox.Text, out double targetMass))
            {
                MessageBox.Show("Неверный формат массы.");
                return;
            }

            reagentPurity.Clear();

            foreach (DataGridViewRow row in dataGridView1.Rows)
            {
                if (row.IsNewRow) continue;

                string compound = row.Cells[0].Value?.ToString();
                string purityStr = row.Cells[1].Value?.ToString();

                if (string.IsNullOrWhiteSpace(compound)) continue;

                if (double.TryParse(purityStr?.Replace(',', '.'), NumberStyles.Any, CultureInfo.InvariantCulture, out double purity))
                {
                    reagentPurity[compound] = purity;
                }
                else
                {
                    reagentPurity[compound] = 100.0; // по умолчанию
                }
            }



            double targetMolarMass = compoundData[targetCompound].MolarMass;
            double targetCoefficient = compoundData[targetCompound].Coefficient;

            double targetMoles = targetMass / targetMolarMass;

            StringBuilder sb = new StringBuilder();
            sb.AppendLine(checkBoxPurity.Checked
                ? "Расчёт с учётом чистоты реагентов:"
                : "Расчёт без учёта чистоты реагентов:");

            foreach (var kvp in compoundData)
            {
                if (kvp.Key == targetCompound)
                    continue;

                double ratio = kvp.Value.Coefficient / targetCoefficient;
                double moles = targetMoles * ratio;
                double mass = moles * kvp.Value.MolarMass;

                if (checkBoxPurity.Checked && reagentPurity.TryGetValue(kvp.Key, out double purity))
                {
                    if (purity <= 0 || purity > 100)
                    {
                        MessageBox.Show($"Неверная чистота для вещества {kvp.Key}. Укажите значение от 0 до 100.");
                        return;
                    }
                    mass /= (purity / 100.0);
                }

                sb.AppendLine($"{kvp.Key}: {mass:F3} г");
            }

            outputMassLog.Text = sb.ToString();
        }

        private void CalculateMasses(string selectedFormula, double knownMass)
        {
            var selected = lastBalancedCompounds.FirstOrDefault(c => c.Formula == selectedFormula);
            if (selected == null)
            {
                MessageBox.Show("Вещество не найдено");
                return;
            }

            double moleCount = knownMass / selected.MolarMass;

            var log = new StringBuilder();
            log.AppendLine($"Из {knownMass:F3} г {selected.Formula} ({moleCount:F4} моль):");

            foreach (var compound in lastBalancedCompounds)
            {
                double ratio = (double)compound.Coefficient / selected.Coefficient;
                double requiredMass = ratio * moleCount * compound.MolarMass;

                log.AppendLine($"{compound.Formula} — {requiredMass:F3} г (коэф. {compound.Coefficient})");
            }

            outputMassLog.Text = log.ToString();
        }

        private void calculateMassButton_Click(object sender, EventArgs e)
        {
            if (productscomboBox.SelectedItem == null || !double.TryParse(massInputBox.Text, out double mass))
            {
                MessageBox.Show("Выберите вещество и введите массу");
                return;
            }

            string formula = productscomboBox.SelectedItem.ToString();

            if (checkBoxPurity.Checked)
            {
                CalculateMassWithPurity(); // учитывает reagentPurity
            }
            else
            {
                CalculateMasses(formula, mass); // обычный расчёт
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            inputBox.Text = "Fe0.8 + O2 + H2 => Fe2O3 + H2O";
            /*
            productscomboBox.Items.Clear();
            foreach (var (compound, i) in products.Select((p, i) => (p, i)))
                productscomboBox.Items.Add($"{compound} (коэфф. {balanceCoefficients[i + reactants.Length]}");
            if (productscomboBox.Items.Count > 0)
                productscomboBox.SelectedIndex = 0;
             */        
            }
    }
}