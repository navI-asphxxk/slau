#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

int gcd(int a, int b){
    for (;;){
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}

int lcm(int a, int b){
    int temp = gcd(a, b);
    return temp ? (a / temp * b) : 0;
}

struct Fraction {
    Fraction() = default;

    Fraction(int num, int denom) : _num(num), _denom(denom) {
        reduceFraction();
    }

    int getNum() const {
        return _num;
    }

    void setNum(int num) {
        _num = num;
        reduceFraction();
    }

    int getDenom() const {
        return _denom;
    }

    void setDenom(int denom) {
        _denom = denom;
        reduceFraction();
    }

    friend std::ostream &operator<<(std::ostream &out, const Fraction &f);

    friend std::istream &operator>>(std::istream &in, Fraction &f);

    Fraction operator+(const Fraction &other) const {
        Fraction res;
        int _lcm = lcm(_denom, other._denom);
        res.setNum(_num * (_lcm / _denom) + other._num * (_lcm / other._denom));
        res.setDenom(_lcm);

        return res;
    }

    Fraction operator+=(const Fraction &other) {
        Fraction res = Fraction{_num, _denom} + other;
        _num = res._num;
        _denom = res._denom;

        return *this;
    }

    Fraction operator-(const Fraction &other) const {
        Fraction res;
        int _lcm = lcm(_denom, other._denom);
        if (!other._denom || !_denom) {

            return {_num, _denom};
        }
        res.setNum(_num * (_lcm / _denom) - other._num * (_lcm / other._denom));
        res.setDenom(_lcm);

        return res;
    }

    Fraction operator*(const Fraction &other) const {
        if (!other._num || !_num) {

            return {0, 1};
        }
        int _gcd1 = gcd(_num, other._denom);
        int _gcd2 = gcd(other._num, _denom);


        int numer1 = _num / _gcd1;
        int denom1 = _denom / _gcd2;

        int numer2 = other._num / _gcd2;
        int denom2 = other._denom / _gcd1;

        return {numer1 * numer2, denom1 * denom2};
    }

    Fraction operator*=(const Fraction &other) {
        Fraction res = {_num * other._num, _denom * other._denom};
        _num = res._num;
        _denom = res._denom;

        return *this;
    }

    Fraction operator/(const Fraction &other) const {

        if (!_num || !other._denom) {
            return {0, 1};
        }


        int numer = _num * other._denom;
        int denom = _denom * other._num;

        int _gcd = gcd(numer, denom);

        return {numer / _gcd, denom / _gcd};

    }

    bool operator>(const Fraction &other) const {
        int _lcm = lcm(_denom, other._denom);

        return _num * (_lcm / _denom) > other._num * (_lcm / other._denom);
    }

    bool operator<(const Fraction &other) const {
        int _lcm = lcm(_denom, other._denom);

        return _num * (_lcm / _denom) < other._num * (_lcm / other._denom);
    }

    bool operator==(const Fraction &other) const {
        int _lcm = lcm(_denom, other._denom);

        return _num * (_lcm / _denom) == other._num * (_lcm / other._denom);
    }

    bool operator!=(const Fraction &other) const {
        int _lcm = lcm(_denom, other._denom);

        return _num * (_lcm / _denom) != other._num * (_lcm / other._denom);
    }

private:
    int _num = 0;
    int _denom = 1;

    void reduceFraction() {
        if (_denom < 0) {
            _num = -_num;
            _denom = -_denom;
        }
        if (_num == 0)
            _denom = 1;

        int _gcd = gcd(_num, _denom);
        _num /= _gcd;
        _denom /= _gcd;
    }
};

std::ostream &operator<<(std::ostream &out, const Fraction &f) {
    if (f._denom == 1)
        out << f._num;
    else
        out << f._num << "/" << f._denom;

    return out;
}

std::istream &operator>>(std::istream &in, Fraction &f) {
    in >> f._num >> f._denom;

    return in;
}

Fraction maxFraction(Fraction f1, Fraction f2) {

    return f1 > f2 ? f1 : f2;
}

Fraction minFraction(Fraction f1, Fraction f2) {

    return f1 < f2 ? f1 : f2;
}

typedef struct equationSystem {
    // хранение матрицы элементов
    std::vector<std::vector<Fraction>> matrix;
    // хранение решений уравнений
    std::vector<Fraction> solutions;

    int nRows;
    int nCols;

    explicit equationSystem() {
        nRows = 0, nCols = 0;
        matrix.resize(nRows);
        solutions.resize(nRows);

        for (int i = 0; i < nRows; i++) {
            matrix[i].resize(nCols);

            for (int j = 0; j < nCols; j++)
                matrix[i][j] = {0, 1};

            solutions[i] = {0,1};
        }
    }

    explicit equationSystem(int rows, int cols) {
        nRows = rows, nCols = cols;
        matrix.resize(nRows);
        solutions.resize(nRows);

        for (int i = 0; i < nRows; i++)
            matrix[i].resize(nCols);

    }

    explicit equationSystem(std::vector<std::vector<Fraction>> &m, std::vector<Fraction>
    &Solutions) {
        nRows = m.size(), nCols = m[0].size();

        matrix.resize(nRows);
        solutions.resize(nRows);

        for (int i = 0; i < nRows; i++)
            matrix[i].resize(nCols);

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++)
                matrix[i][j] = m[i][j];

            solutions[i] = Solutions[i];
        }
    }

    // вывод системы уравнений
    void printSystem() {
        for (int i = 0; i < nRows; i++) {
            bool isFirst = true;
            for (int j = 0; j < nCols; j++)
                if (matrix[i][j] != (Fraction){0,1})
                    continue;
                else if (matrix[i][j] == (Fraction){1,1}) {
                    std::cout << "x" << j + 1;
                    isFirst = false;
                } else if (matrix[i][j] > (Fraction){0,1} && !isFirst) {
                    std::cout << '+' << matrix[i][j] << "x" << j + 1;
                    isFirst = false;
                } else {
                    std::cout << matrix[i][j] << "x" << j + 1;
                    isFirst = false;
                }

            std::cout << " = " << solutions[i] << '\n';
        }

        std::cout << '\n' << std::endl;
    }

    // вывод неупорядоченной системы уравнений
    void printUnorderedSystem(std::vector<int> a) {
        for (int i = 0; i < nRows; i++) {
            bool isFirst = true;

            for (int j = 0; j < nCols; j++)
                if (matrix[i][j] == (Fraction){0,1})
                    continue;
                else if (matrix[i][j] == (Fraction){1,1}) {
                    std::cout << "x" << a[j] + 1;
                    isFirst = false;
                } else if (matrix[i][j] > (Fraction){0,1} && !isFirst) {
                    std::cout << '+' << matrix[i][j] << "x" << a[j] + 1;
                    isFirst = false;
                } else {
                    std::cout << matrix[i][j] << "x" << a[j] + 1;
                    isFirst = false;
                }

            std::cout << " = " << solutions[i] << '\n';
        }

        std::cout << '\n' << std::endl;
    }

    void addRows(int s1, int s2) {
        for (int i = 0; i < nCols; i++)
            matrix[s2][i] += matrix[s1][i];

        solutions[s2] += solutions[s1];
    }

    // Прибавляет к строке индекс row строку b
    void addRow(equationSystem &s, int row, std::vector<Fraction> b) {
        for (int i = 0; i < nCols; i++)
            matrix[row][i] += b[i];
    }

    // Умножает элементы строки row на коэффициент k
    void multipleRows(int row, Fraction k) {
        for (int i = 0; i < nCols; i++)
            matrix[row][i] *= k;
        solutions[row] *= k;
    }

    // Меняет местами строки s1 и s2
    void swapRows(int s1, int s2) {
        std::vector<Fraction> a = matrix[s1];
        matrix[s1] = matrix[s2];
        matrix[s2] = a;
        Fraction b = solutions[s1];
        solutions[s1] = solutions[s2];
        solutions[s2] = b;
    }

    // Возвращает индекс строки с ненулевым элементом после строки с индексом index
    // Если такой строки нет возвращает -1.
    int findNonZeroRow(int index) {
        int numRow = -1;
        int i = index;
        while (numRow < 0 && i < matrix.size()) {
            if (matrix[i][index] != (Fraction){0,1}) {
                numRow = i;
                return numRow;
            } else
                i++;
        }
        return numRow;
    }

    // делает на главной диагонали все элементы ненулевыми
    int makeNonZeroMainDiagonal(equationSystem &source) {
        for (int i = 0; i < source.nRows && source.nCols; i++)
            if (source.matrix[i][i] == (Fraction){0,1}) {
                if (findNonZeroRow(i) != -1) {
                    addRows(findNonZeroRow(i), i);
                    return 1;
                } else
                    return 0;
            }
        return 1;
    }

    // находит и удаляет пустую строку если такая есть
    void deleteEmptyRow(equationSystem &source) {
        for (int i = 0; i < source.nRows; i++) {
            int j = 0;
            while (j < nCols && matrix[i][j] == (Fraction){0,1})
                j++;
            if (j == nCols && solutions[i] == (Fraction){0,1}) {
                nRows--;
                for (int q = i; q < solutions.size() - 1; q++) {
                    matrix[q] = matrix[q + 1];
                    solutions[q] = solutions[q + 1];
                }
                matrix.resize(nRows);
                solutions.resize(nRows);
            }
        }
    }

    Fraction determinant(equationSystem source) {
        int i = 0;
        Fraction determ = matrix[i][i];
        if (determ == (Fraction){0,1})
            return {0,1};
        i++;
        while (i < nRows && i < nCols) {
            if (matrix[i][i] != (Fraction){0,1}) {
                determ *= matrix[i][i];
                i++;
            } else
                return {0,1};
        }
        return determ;
    }

    bool jordanGauss(equationSystem &source) {
        deleteEmptyRow(source);
        if (!makeNonZeroMainDiagonal(source))
            return false;
        for (int i = 0; i < nRows; i++) {
            Fraction mainEl = matrix[i][i];
            Fraction k = (Fraction){1,1} / mainEl;
            multipleRows(i, k);
            for (int j = i + 1; j < nRows; j++) {
                if (matrix[j][i] != (Fraction){0,1}) {
                    Fraction koef = ((Fraction){-1,1} / matrix[j][i]);
                    std::vector<Fraction> strKoef(nCols);
                    for (int qq = 0; qq < nCols; qq++)
                        strKoef[qq] = matrix[i][qq] / koef;
                    addRow(source, j, strKoef);
                    solutions[j] += solutions[i] / koef;
                }
            }
            deleteEmptyRow(source);
            if (!makeNonZeroMainDiagonal(source))
                return false;
        }
        deleteEmptyRow(source);
        if (!makeNonZeroMainDiagonal(source))
            return false;
        for (int i = 1; i < nRows; i++) {
            for (int j = i - 1; j >= 0; j--) {
                if (matrix[j][i] != (Fraction){0,1}) {
                    Fraction koef = (Fraction){-1,1} / matrix[j][i];
                    std::vector<Fraction> strKoef(nCols);
                    for (int k = 0; k < nCols; k++)
                        strKoef[k] = matrix[i][k] / koef;
                    addRow(source, j, strKoef);
                    solutions[j] += solutions[i] / koef;
                }
            }
            deleteEmptyRow(source);
            if (!makeNonZeroMainDiagonal(source))
                return false;
        }
        if (!makeNonZeroMainDiagonal(source)) {
            std::cout << "ERROR";
            return false;
        }
        return true;
    }

    std::vector<Fraction> bazisSolution(equationSystem source, std::vector<int> sequens) {
        std::vector<Fraction> sol(sequens.size(), {0,1});
        for (int i = 0; i < nRows; i++) {
            sol[sequens[i]] = solutions[i];
        }
        return sol;
    }
} SLAU;

bool isRefSolution(std::vector<Fraction> sequens) {
    for (int i = 0; i < sequens.size(); i++)
        if (sequens[i] < (Fraction){0,1})
            return false;
    return true;
}

template<typename Type>
void outputVector(std::vector<Type> a) {
    std::cout << "{";
    for (int i = 0; i < a.size(); i++) {
        std::cout << a[i];
        if (i != a.size() - 1)
            std::cout << "; ";
    }
    std::cout << "}" << std::endl;
}

void inputSLAU1(std::vector<std::vector<Fraction>> &matrix,
                std::vector<Fraction> &Solutions) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++)
            std::cin >> matrix[i][j];
        std::cin >> Solutions[i];
    }
}

//template<typename Type>
void outputSLAU1(std::vector<std::vector<Fraction>> matrix,
                 std::vector<Fraction> Solutions) {
    for (int i = 0; i < matrix.size(); i++) {
        bool isFirst = true;
        for (int j = 0; j < matrix[i].size(); j++)
            if (matrix[i][j] == (Fraction){0,1})
                continue;
            else if (matrix[i][j] == (Fraction){1,1}) {
                std::cout << "x" << j + 1;
                isFirst = false;
            } else if ((matrix[i][j] > (Fraction){1,1}) && !isFirst) {
                std::cout << '+' << matrix[i][j] << "x" << j + 1;
                isFirst = false;
            } else {
                std::cout << matrix[i][j] << "x" << j + 1;
                isFirst = false;
            }
        std::cout << " = " << Solutions[i] << '\n';
    }
}

// Создает сочетания из n по к и записывает их в вектор векторов
void generatingCombinations_(std::vector<int> &a, int k, std::vector<int> &combination, std::vector<bool> &indicator,
                             int sizeCombination, int fix,
                             std::vector<std::vector<int>> &solut) {
    for (int x = fix; x < a.size(); x++) {
        if (indicator.at(x)) {
            combination.at(sizeCombination) = a.at(x);
            indicator.at(x) = false;
            if (sizeCombination + 1 < k) {
                generatingCombinations_(a, k, combination, indicator,
                                        sizeCombination + 1, x, solut);
                indicator.at(x) = true;
            } else {
                solut.push_back(combination);
                indicator.at(x) = true;
            }
        }
    }
}

std::vector<std::vector<int>> generatingCombinations(std::vector<int> &a, int k) {
    static std::vector<int> combination(k);
    static std::vector<bool> indicator(a.size(), true);
    std::vector<std::vector<int>> solut(0, std::vector<int>(0));
    generatingCombinations_(a, k, combination, indicator, 0, 0, solut);

    return solut;
}

bool linearSearch(std::vector<int> a, int b) {
    for (int i = 0; i < a.size(); i++)
        if (a[i] == b)
            return true;
    return false;
}

// Функция x1 + x2 + ... xn
Fraction refFunction(std::vector<Fraction> a) {
    Fraction solution = (Fraction){0,1};
    for (int i = 0; i < a.size(); i++)
        solution += a[i];
    return solution;
}

void inputFraction(Fraction &f, int numer, int denom) {
    f.setNum(numer);
    if (!denom) {
        f.setDenom(1);
    } else {
        f.setDenom(denom);
    }
}

//записывает матрицу matrix и вектор solutions из файла "input.txt",
void inputSLAUFromFile(std::string input, int &r, int &c,
                       std::vector<std::vector<Fraction>> &matrix,
                       std::vector<Fraction> &solutions) {
    std::ifstream inputFile(input);
    std::string s;
    size_t n, m;
    size_t row = 0;
    size_t col = 0;
    bool isFirsIter = true;
    bool isSecondIter = true;

    while (inputFile >> s) {
        if (isFirsIter) {
            n = s[0] - '0';
            isFirsIter = false;
        } else if (!isFirsIter && isSecondIter) {
            m = s[0] - '0';
            matrix.resize(m);
            solutions.resize(n);
            for (size_t q = 0; q < n; q++) {
                matrix[q].resize(n);
            }
            isSecondIter = false;
        } else {
            std::string nsucoor(s);
            std::stringstream ss(nsucoor);
            int numer, denom = 0;
            ss >> numer;
            ss.ignore(1);
            ss >> denom;
            inputFraction(matrix[row][col], numer, denom);
            col++;
            if (col % m == 0) {
                ss >> numer;
                ss.ignore(1);
                ss >> denom;
                inputFraction(solutions[row], numer, denom);
                row++;
                col = 0;
            }
        }
    }
}

int main() {
    std::string s = "input.txt";
    std::ifstream inputFile(s);
    int row = 4;
    int col = 6;

    std::vector<std::vector<Fraction>> matrix(row, std::vector<Fraction>(col));
    std::vector<Fraction> solutions(row);

    matrix = {{{1,1},{-4,1},{8,1},{9,1},{-3,1}, {-1,1}},
              {{8,1},{1,1},{-3,1},{4,1},{5,1}, {6,1}},
              {{4,1},{0,1},{1,1},{3,1},{-2,1}, {-6,1}},
              {{-3,1},{-4,1},{7,1},{6,1},{-1,1}, {4,1}}};
    solutions = {{87,1},{11,1},{17,1},{70,1}};
    equationSystem system1(matrix, solutions);
    equationSystem system2 = system1;
    system2.jordanGauss(system2);
    std::vector<int> a(col, 0);
    for (int i = 0; i < col; i++)
        a[i] = i;

    // генерация всех сочетаний без повторений
    std::vector<std::vector<int>> comb = generatingCombinations(a,
                                                                system2.nRows);
    // перемещения с сочетаниями короче так надо
    std::vector<std::vector<int>> combAll = comb;
    for (int all = 0; all < combAll.size(); all++) {
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < combAll[j].size(); j++) {
                bool check = linearSearch(combAll[all], a[i]);
                if (!check)
                    combAll[all].push_back(a[i]);
            }
        }
    }

    std::vector<std::vector<Fraction>> baseSolutions(0, std::vector<Fraction>(0));
    for (int i = 0; i < combAll.size(); i++) {
        equationSystem currSLAU(system2.nRows, system2.nCols);

        for (int r = 0; r < system2.nRows; r++) {
            for (int c = 0; c < combAll[i].size(); c++)
                currSLAU.matrix[r][c] = system2.matrix[r][combAll[i][c]];
            currSLAU.solutions[r] = system2.solutions[r];
        }
        bool allGood = currSLAU.jordanGauss(currSLAU);
        if ((currSLAU.determinant(currSLAU) != (Fraction){0,1}) && allGood)
            baseSolutions.push_back(currSLAU.bazisSolution(currSLAU, combAll[i]));
    }


    std::vector<std::vector<Fraction>> refSol(0, std::vector<Fraction>(0));
    for (int i = 0; i < baseSolutions.size(); i++)
        if (isRefSolution(baseSolutions[i]))
            refSol.push_back(baseSolutions[i]);
    std::cout << "BASE SOLUTIONS: " << '\n';
    for (int i = 0; i < baseSolutions.size(); i++)
        outputVector(baseSolutions[i]);
    std::cout << '\n' << "REF SOLUTIONS: " << '\n';
    for (int i = 0; i < refSol.size(); i++)
        outputVector(refSol[i]);
    std::string input = " ";


    Fraction maximum = {INT_MIN, 1};
    for (int i = 0; i < baseSolutions.size(); i++)
        maximum = maxFraction(maximum, refFunction(baseSolutions[i]));
    std::cout << "MAX on function is: " << maximum << std::endl;

    Fraction minimum = {INT_MAX, 1};
    for (int i = 0; i < refSol.size(); i++)
        minimum = minFraction(minimum, refFunction(refSol[i]));
    std::cout << "MIN on function is: " << minimum << std::endl;


    return 0;
}