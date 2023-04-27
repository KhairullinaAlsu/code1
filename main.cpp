#include <iostream>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;

const double EPS = 0.00000001;

class ColumnVector{
public:
    int n;
    vector <double> elements;
    ColumnVector(int n) {
        this->n = n;
        elements.resize(n);
    }
    double norm() {
        double res = 0;
        for (int i = 0; i < n; ++i) {
            res += elements[i] * elements[i];
        }
        return sqrt(res);
    }
    ColumnVector operator + (ColumnVector vec) {
        ColumnVector res(this->n);
        for (int i = 0; i < n; i++) {
            res.elements[i] = this->elements[i] + vec.elements[i];
        }
        return res;
    }
    ColumnVector operator - (ColumnVector vec) {
        ColumnVector res(this->n);
        for (int i = 0; i < n; i++) {
            res.elements[i] = this->elements[i] - vec.elements[i];
        }
        return res;
    }
    ColumnVector& operator = (ColumnVector vec) {
        for (int i = 0; i < n; i++) {
            this->elements[i] = vec.elements[i];
        }
        return *this;
    }
};

class Matrix {
public:
    int n, m;
    vector < vector <double> > elements;
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        elements.resize(n, vector <double>(m));
    }
    Matrix operator +(Matrix matrixAdd) {
        Matrix resultAdd(this->n, this->m);
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                resultAdd.elements[i][j] = this->elements[i][j] + matrixAdd.elements[i][j];
            }
        }
        return resultAdd;
    }
    Matrix operator -(Matrix matrixSub) {
        Matrix resultSub(this->n, this->m);
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                resultSub.elements[i][j] = this->elements[i][j] - matrixSub.elements[i][j];
            }
        }
        return resultSub;
    }
    Matrix operator *(Matrix matrixMul) {
        Matrix resultMul(this->n,  matrixMul.m);
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < matrixMul.m; j++) {
                for (int l = 0; l < min(this->m, matrixMul.n); l++) {
                    resultMul.elements[i][j] += this->elements[i][l] * matrixMul.elements[l][j];
                }
            }
        }
        return resultMul;
    }
    ColumnVector operator * (ColumnVector vec) {
        ColumnVector res(this->n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res.elements[i] += this->elements[i][j] * vec.elements[j];
            }
        }
        return res;
    }
    virtual Matrix& operator =(Matrix a) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; j++) {
                elements[i][j] = a.elements[i][j];
            }
        }
        return *this;
    }
    Matrix transpose() {
        Matrix resultTrans(this->m, this->n);
        for (int i = 0; i < this->m; i++) {
            for (int j = 0; j < this->n; j++) {
                resultTrans.elements[i][j] = this->elements[j][i];
            }
        }
        return resultTrans;
    }
};

istream& operator >> (istream& in, Matrix& matrix) {
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.m; j++) {
            in >> matrix.elements[i][j];
        }
    }
    return in;
}

ostream& operator << (ostream& out, Matrix& matrix) {
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.m; j++) {
            if (abs(matrix.elements[i][j]) < EPS) {
                out << 0.0 << " ";
            }
            else {
                out << matrix.elements[i][j] << " ";
            }
        }
        out << '\n';
    }
    return out;
}

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) { }
    SquareMatrix& operator =(Matrix matrixEq) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                this->elements[i][j] = matrixEq.elements[i][j];
            }
        }
        return *this;
    }
};

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            elements[i][i] = 1;
        }
    }
};

istream& operator >> (istream& in, ColumnVector& matrix) {
    for (int i = 0; i < matrix.n; i++) {
        in >> matrix.elements[i];
    }
    return in;
}

ostream& operator << (ostream& out, ColumnVector& matrix) {
    for (int i = 0; i < matrix.n; i++) {
        out << fixed << setprecision(4) << matrix.elements[i] << '\n';
    }
    return out;
}

class PermutationMatrix : public Matrix {
public:
    PermutationMatrix(Matrix* A, int n, int i, int j) : Matrix(n, n) {
        for (int l = 0; l < n; l++) {
            elements[l][l] = 1;
        }
        for (int l = 0; l < n; l++) {
            double tmp = elements[i][l];
            elements[i][l] = elements[j][l];
            elements[j][l] = tmp;
        }
    }
};

class EliminationMatrix : public Matrix {
public:
    EliminationMatrix(Matrix* A, int n, int i, int j) : Matrix(A->n, A->m) {
        double e = -1.0 * (A->elements[i - 1][j - 1] / A->elements[j - 1][j - 1]);
        for (int l = 0; l < A->n; l++) {
            for (int k = 0; k < A->m; k++) {
                if (l == k) {
                    this->elements[l][k] = 1;
                }
                else if (l == i - 1 and k == j - 1) {
                    this->elements[l][k] = e;
                }
                else {
                    this->elements[l][k] = 0;
                }
            }
        }
    }
};

class AugmentedMatrix : public Matrix {
public:
    AugmentedMatrix(Matrix* a, int n) : Matrix(a->n, a->n * 2) {
        for (int i = 0; i < a->n; ++i) {
            for (int j = 0; j < a->n; ++j) {
                this->elements[i][j] = a->elements[i][j];
            }
            this->elements[i][i + n] = 1;
        }
    }
};

void DirectWay(Matrix& a) {
    for (int i = 0; i < a.n; i++) {
        int max_row = i;
        for (int j = i; j < a.n; j++) {
            if (EPS + abs(a.elements[j][i]) > abs(a.elements[max_row][i]) && a.elements[j][i] != a.elements[max_row][i]) {
                max_row = j;
            }
        }
        if(abs(a.elements[max_row][i]) < EPS) continue;
        if (max_row != i) {
            PermutationMatrix permutationMatrix(&a, a.n, i, max_row);
            a = permutationMatrix * a;
        }
        for (int j = i + 1; j < a.n; j++) {
            if (abs(a.elements[j][i]) < EPS) {
                continue;
            }
            EliminationMatrix eliminationMatrix(&a, a.n, j + 1, i + 1);
            a = eliminationMatrix * a;
        }
    }
}

void WayBack(Matrix& a) {
    for (int i = a.n - 1; i >= 0; i--) {
        double factor = a.elements[i][i];
        for (int j = i; j < a.n * 2; j++) {
            if (abs(factor) < EPS){
                continue;
            }
            a.elements[i][j] /= factor;
        }
        for (int j = i - 1; j >= 0; j--) {
            if (abs(a.elements[j][i]) < EPS) {
                continue;
            }
            EliminationMatrix eliminationMatrix(&a, a.n, j + 1, i + 1);
            a = eliminationMatrix * a;
        }
    }
}

Matrix isolation(Matrix& a) {
    Matrix diagonalMatrixInverse(a.n, a.n);
    for (int i = 0; i < a.n; i++) {
        for (int j = 0; j < a.n; j++) {
            diagonalMatrixInverse.elements[i][j] = a.elements[i][j + a.n];
        }
    }
    return diagonalMatrixInverse;
}


Matrix matrixInverse(Matrix& a) {
    AugmentedMatrix augmentedMatrix(&a, a.n);
    DirectWay(augmentedMatrix);
    WayBack(augmentedMatrix);
    //normalization(augmentedMatrix);
    return isolation(augmentedMatrix);
}

int main() {
    cout << fixed << setprecision(4);
    int m, n;
    cin >> m;
    ColumnVector t(m);
    ColumnVector b(m);
    for (int i = 0; i < m; ++i) {
        cin >> t.elements[i] >> b.elements[i];
    }
    cin >> n;
    Matrix A(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            A.elements[i][j] = pow(t.elements[i], j);
        }
    }
    Matrix ATranspose = A.transpose();
    Matrix ATransposeMulA = ATranspose * A;
    Matrix ATransposeMulAInv = matrixInverse(ATransposeMulA);
    ColumnVector ATransposeMulB = ATranspose * b;
    ColumnVector x = ATransposeMulAInv * ATransposeMulB;
    cout << "A:\n" << A << "A_T*A:\n" << ATransposeMulA << "(A_T*A)^-1:\n" << ATransposeMulAInv << "A_T*b:\n" << ATransposeMulB << "x~:\n" << x;
    return 0;
}