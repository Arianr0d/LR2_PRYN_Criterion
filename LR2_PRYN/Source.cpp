#include <iostream>;
#include <sstream>;
#include <fstream>;
#include <vector>;
#include <algorithm>;
#include <string>;
#include <cmath>;

using namespace std;

/*
	Поиск решения игры в условиях неопределенности, применяя критерии Вальда, Сэвиджа, Гурвица,
	Байеса-Лапласа и Ходжа-Лемана
*/

// функция записи результата в файл
void output(vector<pair<int, int>> res_wald, vector<pair<int, int>> res_savage, vector<pair<int, float>> res_hurwitz,
	vector<pair<int, float>> res_bayeslaplace, vector<pair<int, float>> res_hodgeleman) {

	ofstream out;
	out.open("output.txt");
	if (out.is_open()) {
		out << "решения игры следующие:" << endl;
		out << "по критерию Вальда: ";
		for (int i = 0; i < res_wald.size(); i++) {
			out << to_string(res_wald[i].second) + " при стратегии " +
				to_string(res_wald[i].first) << endl;
		}
		out << "по критерию Сэвиджа: ";
		for (int i = 0; i < res_savage.size(); i++) {
			out << to_string(res_savage[i].second) + " при стратегии " +
				to_string(res_savage[i].first) << endl;
		}
		out << "по критерию Гурвица: ";
		for (int i = 0; i < res_hurwitz.size(); i++) {
			out << res_hurwitz[i].second;
			out << " при стратегии ";
			out << to_string(res_hurwitz[i].first) << endl;
		}
		out << "по критерию Байеса–Лапласа: ";
		for (int i = 0; i < res_bayeslaplace.size(); i++) {
			out << res_bayeslaplace[i].second;
			out << " при стратегии ";
			out << res_bayeslaplace[i].first << endl;
		}
		out << "по критерию Ходжа–Лемана: ";
		for (int i = 0; i < res_hodgeleman.size(); i++) {
			out << res_hodgeleman[i].second;
			out << " при стратегии ";
			out << res_hodgeleman[i].first << endl;
		}
	}
}

// математического ожидания
float findMathExpectation(vector<int> &row) {

	float sum = 0;
	for (int i = 0; i < row.size(); i++) {
		sum += row[i];
	}
	sum /= row.size();
	return sum;
}

// функция поиска матрицы рисков
void findMatrixRisk(vector<int> max_column, vector<vector<int>>& matrix_risk) {

	for (int i = 0; i < matrix_risk[0].size(); i++) {
		for (int j = 0; j < matrix_risk.size(); j++) {
			matrix_risk[j][i] = max_column[i] - matrix_risk[j][i];
		}
	}
}

// функция поиска всех максимумов по вектору
template <typename T1, typename T2>
void findAllMaximum(vector<T1> vec, vector<pair<int, T2>> &result, bool max_value = true) {

	/*auto it = minmax_element(vec.begin(), vec.end());
	int min_idx = distance(vec.begin(), it.first);
	int max_idx = distance(vec.begin(), it.second);*/

	/*
		Выводить один макс / мин или все найденные!?
	*/

	T2 minmax_value = vec[0];
	result.push_back(make_pair(1, minmax_value));

	for (int i = 1; i < vec.size(); i++) {
		if (minmax_value < vec[i] && max_value) {
			minmax_value = vec[i];
			result.clear();
			result.push_back(make_pair(i + 1, minmax_value));
		}
		else if (minmax_value > vec[i] && !max_value) {
			minmax_value = vec[i];
			result.clear();
			result.push_back(make_pair(i + 1, minmax_value));
		}
		else if (minmax_value == vec[i]) {
			result.push_back(make_pair(i + 1, minmax_value));
		}
	}
}

// функция поиска минимума по вектору (строка / столбец)
template <typename T>
int findMinMax(vector<T> row, bool max_value = true) {

	if (row.size() == 0) {
		cout << "Пустой массив!" << endl;
		return 0;
	}
	
	sort(row.begin(), row.end());
	if (!max_value) return row.front();
	else return row.back();
}

// функция поиска оптимальной стратегии по критерию Ходжа-Лемана
void findHodge_LemanCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result, float mu = 0.5) {

	vector<float> W;
	// поиск среднего арифметического минимума строки и ее математического ожидания
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMathExpectation(matrix[i]) * mu + findMinMax(matrix[i], false) * (1 - mu));
	}

	// поиск максимальной стратегии
	findAllMaximum(W, result, true);
}

// функция поиска оптимальной стратегии по критерию Байеса-Лапласа
void findBayes_LaplaceCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result) {

	vector<float> math_exp;
	// расчёт математического ожидания каждой строки
	for (int i = 0; i < matrix.size(); i++) {
		math_exp.push_back(findMathExpectation(matrix[i]));
	}

	// поиск максимальной стратегии
	findAllMaximum(math_exp, result, true);
}

// функция поиска оптимальной стратегии по критерию Гурвица
void findHurwitzCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result, float lambda = 0.5) {

	vector<float> W;
	// вычисление значения критерия
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMinMax(matrix[i], true) * lambda + findMinMax(matrix[i], false) * (1 - lambda));
	}

	// поиск максимальной стратегии
	findAllMaximum(W, result, true);
}

// функция поиска оптимальной стратегии по критерию Сэвиджа
void findSavageCriterion(vector<vector<int>>& matrix, vector<pair<int, int>>& result) {

	vector<int> D, D_column;
	vector<vector<int>> m_risk = matrix;

	// поиск максимального значения по столбцам
	for (int i = 0; i < matrix[0].size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			D_column.push_back(matrix[j][i]);
		}
		D.push_back(findMinMax(D_column, true));
		D_column.clear();
	}

	// поиск матрицы рисков R
	findMatrixRisk(D, std::ref(m_risk));
	D.clear();
	// поиск максимума по строкам матрицы рисков R
	for (int i = 0; i < m_risk.size(); i++) {
		D.push_back(findMinMax(m_risk[i], true));
	}
	// поиск минимума среди максимумов
	findAllMaximum(D, result, false);
}

// функция поиска оптимальной стратегии по критерию Вальда
void findWaldCriterion(vector<vector<int>> &matrix, vector<pair<int, int>> &result) {

	vector<int> W;
	// поиск минимального значения по строке
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMinMax(matrix[i], false));
	}

	// поиск максимума среди минимумов
	findAllMaximum(W, result, true);
}

// функция считывания матрицы из консоли
void inputMatrix(int n, int m, vector<vector<int>> &matrix) {

	cin.ignore(cin.rdbuf()->in_avail());
	for (int i = 0; i < n; i++) {
		int j = 0;
		string row;
		getline(cin, row);
		istringstream ist(row);

		while (ist >> row) {
			matrix[i][j] = stoi(row);
			j++;
		}
	}
}

// функция считывания матрицы из файла
void importFromFile(int n, int m, vector<vector<int>> &matrix) {

	string str;
	ifstream in("input.txt");
	if (in.is_open()) {
		int i = 0;
		while (getline(in, str)) {

			if (str.size() == 0) continue;

			string row;
			istringstream ist(str);
			int j = 0;

			while (ist >> row) {
				matrix[i][j] = stoi(row);
				j++;
			}
			i++;
		}
	}
	in.close();
}

int main() {

	setlocale(LC_ALL, "ru");

	string ans;
	float lambda = 0.5, mu = 0.5;
	int n = 8, m = 10;
	vector<vector<int>> A(n, vector<int>(m));

	cout << "Поиск решения игры в условиях неопределенности, применяя критерии Вальда, Сэвиджа,\n" <<
		"Гурвица, Байеса - Лапласа и Ходжа - Лемана" << endl << endl;
	while (true) {
		cout << "Использовать значения по умолчанию (2 вариант): (yes / no) ";
		cin >> ans;
		if (ans == "no" || ans == "yes") break;
	}

	if (ans == "no") {
		// считывание данных с консоли
		cout << "Введите число строк матрицы A: ";
		cin >> n;
		cout << "Введите число столбцов матрицы A: ";
		cin >> m;
		cout << "Введите матрицу A: " << endl;
		inputMatrix(n, m, std::ref(A));
	}
	else {
		// считывание данных из файла
		importFromFile(n, m, std::ref(A));
	}

	vector<pair<int, int>> res_wald;
	// поиск оптимальной стратегии по критерию Вальда
	findWaldCriterion(A, std::ref(res_wald));

	vector<pair<int, int>> res_savage;
	// поиск оптимальной стратегии по критерию Сэвиджа
	findSavageCriterion(A, std::ref(res_savage));

	vector<pair<int, float>> res_hurwitz;
	// поиск оптимальной стратегии по критерию Гурвица
	findHurwitzCriterion(A, std::ref(res_hurwitz), lambda);

	vector<pair<int, float>> res_bayeslaplace;
	// поиск оптимальной стратегии по критерию Байеса-Лапласа
	findBayes_LaplaceCriterion(A, std::ref(res_bayeslaplace));

	vector<pair<int, float>> res_hodgeleman;
	// поиск оптимальной стратегии по критерию Ходжа-Лемана
	findHodge_LemanCriterion(A, std::ref(res_hodgeleman), mu);

	// запись результата в файл
	output(res_wald, res_savage, res_hurwitz, res_bayeslaplace, res_hodgeleman);

	cout << "Результаты записаны в файл!";
}