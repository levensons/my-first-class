from sys import stdin
import copy


def fast_pow(my_matrix, k):
    n, m = my_matrix.size()
    if k == 0:
        res_matrix = [] * n
        for i in range(n):
            res_matrix.append([])
        for i in range(n):
            for j in range(m):
                if i == j:
                    res_matrix[i].append(1)
                else:
                    res_matrix[i].append(0)
        return Matrix(res_matrix)
    if k % 2 == 1:
        return my_matrix * fast_pow(my_matrix, k - 1)
    res = fast_pow(my_matrix, k // 2)
    return res * res


def NOD(a, b):
    if not a or not b:
        return a + b
    return NOD(b, a % b)


class MatrixError(BaseException):
    def __init__(self, matrix, other):
        self.matrix1 = matrix
        self.matrix2 = other


class MatrixError1(BaseException):
    def __init__(self):
        self.matrix = [[1]]


class Matrix:
    def __init__(self, matrix):
        new_matrix = copy.deepcopy(matrix)
        self.matrix = new_matrix

    def __str__(self):
        my_str = ''
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                my_str += str(self.matrix[i][j])
                if j != len(self.matrix[i]) - 1:
                    my_str += '\t'
            if i != len(self.matrix) - 1:
                my_str += '\n'
        return my_str

    def size(self):
        return len(self.matrix), len(self.matrix[0])

    def __add__(self, other):
        new_matrix = []
        n, m = other.size()
        n1, m1 = self.size()
        if n != n1 or m != m1:
            raise MatrixError(self, other)
        for i in range(n):
            new_matrix.append([])
        for i in range(n):
            for j in range(m):
                new_matrix[i].append(self.matrix[i][j] + other.matrix[i][j])
        return Matrix(new_matrix)

    def __mul__(self, other):
        n1, m1 = self.size()
        new_matrix = []
        for i in range(n1):
            new_matrix.append([])
        if isinstance(other, int) or isinstance(other, float):
            for i in range(n1):
                for j in range(m1):
                    new_matrix[i].append(self.matrix[i][j] * other)
        elif isinstance(other, Matrix):
            n2, m2 = other.size()
            if m1 != n2:
                raise MatrixError(self, other)
            for i in range(n1):
                for j in range(m2):
                    res = 0
                    for k in range(m1):
                        res += self.matrix[i][k] * other.matrix[k][j]
                    new_matrix[i].append(res)
        else:
            raise MatrixError(self, other)
        return Matrix(new_matrix)

    __rmul__ = __mul__

    def transpose(self):
        n, m = self.size()
        new_matrix = []
        for i in range(m):
            new_matrix.append([])
        for i in range(n):
            for j in range(m):
                new_matrix[j].append(self.matrix[i][j])
        self.matrix = copy.deepcopy(new_matrix)
        return self

    def solve(self, vect_arr):
        n, m = self.size()
        if m > n:
            raise MatrixError(self, [1])
        m += 1
        new_matrix = []
        for i in range(n):
            new_matrix.append([])
        for i in range(n):
            for j in range(m - 1):
                new_matrix[i].append((self.matrix[i][j], 1))
        for i in range(n):
            new_matrix[i].append((vect_arr[i], 1))
        # привожу к ступенчатому виду
        for i in range(n - 1):
            j = 0
            while new_matrix[i][j] == (0, 0) and j < m - 1:
                j += 1
            if j == m - 1 and i != n - 1:
                break
            a0 = new_matrix[i][j]
            for k in range(i + 1, n - 1):
                ak = new_matrix[k][j]
                new_matrix[k][j] = (0, 0)
                for t in range(j + 1, m):
                    at = new_matrix[k][t]
                    ai = new_matrix[i][t]
                    p1 = ai[1] * at[1] * a0[0] * ak[1]
                    p2 = ai[0] * a0[1] * ak[0] * at[1]
                    p3 = at[0] * a0[0] * ak[1] * ai[1]
                    new_matrix[k][t] = (p3 - p2, p1)
                    Nod = NOD(new_matrix[k][t][0], new_matrix[k][t][1])
                    if Nod:
                        a = new_matrix[k][t][0] // Nod
                        b = new_matrix[k][t][1] // Nod
                        new_matrix[k][t] = (a, b)
                    if new_matrix[k][t][0] == 0:
                        new_matrix[k][t] = (0, 0)
        # превращаю матрицу в обычный вид
        for i in range(n):
            for j in range(m):
                if new_matrix[i][j] != (0, 0):
                    a = new_matrix[i][j][0]
                    b = new_matrix[i][j][1]
                    new_matrix[i][j] = a / b
                else:
                    new_matrix[i][j] = 0
        # проверяю на нулевую строку с ненулевый свободным членом
        for i in range(n):
            flag = 1
            for j in range(m - 1):
                if new_matrix[i][j]:
                    flag = 0
                    break
            if flag and new_matrix[i][m - 1]:
                raise MatrixError1(self)
        # ищу нулевые строки
        zero = [0] * n
        for i in range(n):
            cnt = 0
            for j in range(m):
                if new_matrix[i][j] == 0:
                    cnt += 1
            if cnt == m:
                zero[i] = 1
        # ищу главные переменные
        main_columns = [0] * m
        col_idx = 0
        for i in range(n):
            while new_matrix[i][col_idx] == 0 and col_idx < m - 1:
                col_idx += 1
            main_columns[col_idx] = i
        # проверяю что во всех столбцах есть главная переменная
        # это значит что кол-во ур-ий равно кол-ву неизв.
        sum1 = 0
        sum2 = 0
        for i in range(m - 1):
            sum1 += main_columns[i]
        for i in range(n):
            sum2 += zero[i]
        sum2 = n - sum2
        if sum1 != sum2:
            raise MatrixError1(self)
        # получаю ответ
        col = m - 2
        answer = []
        for i in range(n - 1, -1, -1):
            if not zero[i]:
                res = new_matrix[i][m - 1]
                for j in range(col + 1, m - 1):
                    res -= new_matrix[i][j] * answer[m - j - 2]
                answer.append(res / new_matrix[i][col])
                col -= 1
        answer.reverse()
        return answer

    @staticmethod
    def transposed(matrix):
        n, m = matrix.size()[0], matrix.size()[1]
        new_matrix = []
        for i in range(m):
            new_matrix.append([])
        for i in range(n):
            for j in range(m):
                new_matrix[j].append(matrix.matrix[i][j])
        return Matrix(new_matrix)


class SquareMatrix(Matrix):
    def __pow__(self, k):
        return fast_pow(self, k)

exec(stdin.read())
