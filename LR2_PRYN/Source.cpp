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
		if (!res_wald.empty()) out << to_string(res_wald[0].second) + " ��� ��������� " +
				to_string(res_wald[0].first) << endl;
		out << "�� �������� �������: ";
		if (!res_savage.empty()) out << to_string(res_savage[0].second) + " ��� ��������� " +
				to_string(res_savage[0].first) << endl;
		out << "�� �������� �������: ";
		if (!res_hurwitz.empty()) {
			out << res_hurwitz[0].second;
			out << " ��� ��������� ";
			out << to_string(res_hurwitz[0].first) << endl;
		}
		out << "�� �������� ��������������: ";
		if (!res_bayeslaplace.empty()) {
			out << res_bayeslaplace[0].second;
			out << " ��� ��������� ";
			out << res_bayeslaplace[0].first << endl;
		}
		out << "�� �������� ������������: ";
		if (!res_hodgeleman.empty()) {
			out << res_hodgeleman[0].second;
			out << " ��� ��������� ";
			out << res_hodgeleman[0].first << endl;
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

// ������� ���������� ������� �� �����
void importFromFile(int &n, int &m, float &lambda, float &mu, vector<vector<int>> &matrix) {

	string str;
	ifstream in("input.txt");
	if (in.is_open()) {
		int i = 0, k = 0;
		while (getline(in, str)) {

			if (str.size() == 0) continue;
			
			switch (k) {
				case 0: 
					n = stoi(str); // ����� ����� � �������
					break;
				case 1: 
					m = stoi(str);  // ����� �������� � �������
					break;
				case 2: 
					lambda = stof(str); // �������� ������ ��� �������� ������� 
					break;
				case 3: 
					mu = stof(str); // �������� �� ��� �������� �����-������
					break;

				default: 
					string row;
					istringstream ist(str);
					int j = 0;

					while (ist >> row) {
						matrix[i][j] = stoi(row);
						j++;
					}
					i++;
					break;
			}
			k++;
		}
	}
	in.close();
}

int main() {

	setlocale(LC_ALL, "ru");

	float lambda = 0.5, mu = 0.5;
	int n = 8, m = 10;
	vector<vector<int>> A(n, vector<int>(m));

	cout << "����� ������� ���� � �������� ����������������, �������� �������� ������, �������,\n" <<
		"�������, ������ - ������� � ����� - ������" << endl << endl;
	cout << "���������� ������ �� �����!" << endl;

	// ���������� ������ �� �����
	importFromFile(n, m, lambda, mu, ref(A));

	vector<pair<int, int>> res_wald;
	// ����� ����������� ��������� �� �������� ������
	findWaldCriterion(A, ref(res_wald));

	vector<pair<int, int>> res_savage;
	// ����� ����������� ��������� �� �������� �������
	findSavageCriterion(A, ref(res_savage));

	vector<pair<int, float>> res_hurwitz;
	// ����� ����������� ��������� �� �������� �������
	findHurwitzCriterion(A, ref(res_hurwitz), lambda);

	vector<pair<int, float>> res_bayeslaplace;
	// ����� ����������� ��������� �� �������� ������-�������
	findBayes_LaplaceCriterion(A, ref(res_bayeslaplace));

	vector<pair<int, float>> res_hodgeleman;
	// ����� ����������� ��������� �� �������� �����-������
	findHodge_LemanCriterion(A, ref(res_hodgeleman), mu);

	// ������ ���������� � ����
	output(res_wald, res_savage, res_hurwitz, res_bayeslaplace, res_hodgeleman);

	cout << "���������� ������� ���������!" << endl;
}