#include <iostream>;
#include <sstream>;
#include <fstream>;
#include <vector>;
#include <algorithm>;
#include <string>;
#include <cmath>;

using namespace std;

/*
	����� ������� ���� � �������� ����������������, �������� �������� ������, �������, �������,
	������-������� � �����-������
*/

// ������� ������ ���������� � ����
void output(vector<pair<int, int>> res_wald, vector<pair<int, int>> res_savage, vector<pair<int, float>> res_hurwitz,
	vector<pair<int, float>> res_bayeslaplace, vector<pair<int, float>> res_hodgeleman) {

	ofstream out;
	out.open("output.txt");
	if (out.is_open()) {
		out << "������� ���� ���������:" << endl;
		out << "�� �������� ������: ";
		for (int i = 0; i < res_wald.size(); i++) {
			out << to_string(res_wald[i].second) + " ��� ��������� " +
				to_string(res_wald[i].first) << endl;
		}
		out << "�� �������� �������: ";
		for (int i = 0; i < res_savage.size(); i++) {
			out << to_string(res_savage[i].second) + " ��� ��������� " +
				to_string(res_savage[i].first) << endl;
		}
		out << "�� �������� �������: ";
		for (int i = 0; i < res_hurwitz.size(); i++) {
			out << res_hurwitz[i].second;
			out << " ��� ��������� ";
			out << to_string(res_hurwitz[i].first) << endl;
		}
		out << "�� �������� ��������������: ";
		for (int i = 0; i < res_bayeslaplace.size(); i++) {
			out << res_bayeslaplace[i].second;
			out << " ��� ��������� ";
			out << res_bayeslaplace[i].first << endl;
		}
		out << "�� �������� ������������: ";
		for (int i = 0; i < res_hodgeleman.size(); i++) {
			out << res_hodgeleman[i].second;
			out << " ��� ��������� ";
			out << res_hodgeleman[i].first << endl;
		}
	}
}

// ��������������� ��������
float findMathExpectation(vector<int> &row) {

	float sum = 0;
	for (int i = 0; i < row.size(); i++) {
		sum += row[i];
	}
	sum /= row.size();
	return sum;
}

// ������� ������ ������� ������
void findMatrixRisk(vector<int> max_column, vector<vector<int>>& matrix_risk) {

	for (int i = 0; i < matrix_risk[0].size(); i++) {
		for (int j = 0; j < matrix_risk.size(); j++) {
			matrix_risk[j][i] = max_column[i] - matrix_risk[j][i];
		}
	}
}

// ������� ������ ���� ���������� �� �������
template <typename T1, typename T2>
void findAllMaximum(vector<T1> vec, vector<pair<int, T2>> &result, bool max_value = true) {

	/*auto it = minmax_element(vec.begin(), vec.end());
	int min_idx = distance(vec.begin(), it.first);
	int max_idx = distance(vec.begin(), it.second);*/

	/*
		�������� ���� ���� / ��� ��� ��� ���������!?
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

// ������� ������ �������� �� ������� (������ / �������)
template <typename T>
int findMinMax(vector<T> row, bool max_value = true) {

	if (row.size() == 0) {
		cout << "������ ������!" << endl;
		return 0;
	}
	
	sort(row.begin(), row.end());
	if (!max_value) return row.front();
	else return row.back();
}

// ������� ������ ����������� ��������� �� �������� �����-������
void findHodge_LemanCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result, float mu = 0.5) {

	vector<float> W;
	// ����� �������� ��������������� �������� ������ � �� ��������������� ��������
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMathExpectation(matrix[i]) * mu + findMinMax(matrix[i], false) * (1 - mu));
	}

	// ����� ������������ ���������
	findAllMaximum(W, result, true);
}

// ������� ������ ����������� ��������� �� �������� ������-�������
void findBayes_LaplaceCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result) {

	vector<float> math_exp;
	// ������ ��������������� �������� ������ ������
	for (int i = 0; i < matrix.size(); i++) {
		math_exp.push_back(findMathExpectation(matrix[i]));
	}

	// ����� ������������ ���������
	findAllMaximum(math_exp, result, true);
}

// ������� ������ ����������� ��������� �� �������� �������
void findHurwitzCriterion(vector<vector<int>>& matrix, vector<pair<int, float>>& result, float lambda = 0.5) {

	vector<float> W;
	// ���������� �������� ��������
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMinMax(matrix[i], true) * lambda + findMinMax(matrix[i], false) * (1 - lambda));
	}

	// ����� ������������ ���������
	findAllMaximum(W, result, true);
}

// ������� ������ ����������� ��������� �� �������� �������
void findSavageCriterion(vector<vector<int>>& matrix, vector<pair<int, int>>& result) {

	vector<int> D, D_column;
	vector<vector<int>> m_risk = matrix;

	// ����� ������������� �������� �� ��������
	for (int i = 0; i < matrix[0].size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			D_column.push_back(matrix[j][i]);
		}
		D.push_back(findMinMax(D_column, true));
		D_column.clear();
	}

	// ����� ������� ������ R
	findMatrixRisk(D, std::ref(m_risk));
	D.clear();
	// ����� ��������� �� ������� ������� ������ R
	for (int i = 0; i < m_risk.size(); i++) {
		D.push_back(findMinMax(m_risk[i], true));
	}
	// ����� �������� ����� ����������
	findAllMaximum(D, result, false);
}

// ������� ������ ����������� ��������� �� �������� ������
void findWaldCriterion(vector<vector<int>> &matrix, vector<pair<int, int>> &result) {

	vector<int> W;
	// ����� ������������ �������� �� ������
	for (int i = 0; i < matrix.size(); i++) {
		W.push_back(findMinMax(matrix[i], false));
	}

	// ����� ��������� ����� ���������
	findAllMaximum(W, result, true);
}

// ������� ���������� ������� �� �������
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

// ������� ���������� ������� �� �����
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

	cout << "����� ������� ���� � �������� ����������������, �������� �������� ������, �������,\n" <<
		"�������, ������ - ������� � ����� - ������" << endl << endl;
	while (true) {
		cout << "������������ �������� �� ��������� (2 �������): (yes / no) ";
		cin >> ans;
		if (ans == "no" || ans == "yes") break;
	}

	if (ans == "no") {
		// ���������� ������ � �������
		cout << "������� ����� ����� ������� A: ";
		cin >> n;
		cout << "������� ����� �������� ������� A: ";
		cin >> m;
		cout << "������� ������� A: " << endl;
		inputMatrix(n, m, std::ref(A));
	}
	else {
		// ���������� ������ �� �����
		importFromFile(n, m, std::ref(A));
	}

	vector<pair<int, int>> res_wald;
	// ����� ����������� ��������� �� �������� ������
	findWaldCriterion(A, std::ref(res_wald));

	vector<pair<int, int>> res_savage;
	// ����� ����������� ��������� �� �������� �������
	findSavageCriterion(A, std::ref(res_savage));

	vector<pair<int, float>> res_hurwitz;
	// ����� ����������� ��������� �� �������� �������
	findHurwitzCriterion(A, std::ref(res_hurwitz), lambda);

	vector<pair<int, float>> res_bayeslaplace;
	// ����� ����������� ��������� �� �������� ������-�������
	findBayes_LaplaceCriterion(A, std::ref(res_bayeslaplace));

	vector<pair<int, float>> res_hodgeleman;
	// ����� ����������� ��������� �� �������� �����-������
	findHodge_LemanCriterion(A, std::ref(res_hodgeleman), mu);

	// ������ ���������� � ����
	output(res_wald, res_savage, res_hurwitz, res_bayeslaplace, res_hodgeleman);

	cout << "���������� �������� � ����!";
}