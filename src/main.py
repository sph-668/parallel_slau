import time

class SparseMatrix:
    
    def __init__(self):
        self.n = 0
        self.ai = []
        self.aj = []
        self.di = []
        self.al = []
        self.au = []
    
    def input(self, file):
        self.n = int(file.readline())
        if (self.n>0):
            self.ai=list(map(int,file.readline().split()))
            if (self.ai[self.n-1] > 0):
                self.aj=list(map(int,file.readline().split()))
            self.di = list(map(float,file.readline().split()))
            if (self.ai[self.n-1] > 0):
                self.al = list(map(float, file.readline().split()))
                self.au = list(map(float, file.readline().split()))
    
    def factorizeLU(self):
        for i in range(0,self.n):
            sd = 0.
            i0 = self.ai[i]
            i1 = self.ai[i+1]
            for j in range(i0,i1):
                sl = 0.
                su = 0.
                j0 = self.ai[self.aj[j]]
                j1 = self.ai[self.aj[j]+1]
                ki = i0
                kj = j0
                while (ki < j and kj < j1):
                    jl = self.aj[ki]
                    ju = self.aj[kj]
                    if (jl == ju):
                        sl += self.au[kj] * self.al[ki]
                        su += self.al[kj] * self.au[ki]
                        ki+=1
                        kj+=1
                    elif (jl < ju):
                        ki+=1
                    else:
                        kj+=1
                self.au[j] = (self.au[j] - su) / self.di[self.aj[j]]
                self.al[j] = (self.al[j] - sl)
                sd += self.au[j] * self.al[j]
            self.di[i] = self.di[i] - sd
            if (abs(self.di[i] <= 1e-15)):
                return False
        return True
    
    def get(self, i, j):
        if (i==j):
            return True, self.di[i]
        v = self.ai if i>j else self.au
        if (i < j):
            i, j = j, i
        _i = self.ai[i]
        while (_i < self.ai[i+1] and self.aj[_i]<j):
            _i+=1
        if (_i < self.ai[i+1] and j == self.aj[_i]):
            return True, v[_i]
        return False, 0.
    
    def set(self, i, j, val):
        if (i==j):
            self.di[i] = val
            return True
        v = self.ai if i>j else self.au
        if (i < j):
            i, j = j, i
        _i = self.ai[i]
        while (_i < self.ai[i+1] and self.aj[_i]<j):
            _i+=1
        if (_i < self.ai[i+1] and j == self.aj[_i]):
            v[_i] = val
            return True
        return False


def find_y(sm, b):
        for i in range(0, sm.n):
            sum = 0
            for idx in range(sm.ai[i],sm.ai[i+1]):
                sum += sm.al[idx] * b[sm.aj[idx]]
            b[i] = (b[i] - sum) / sm.di[i]
    
def find_x(sm, y):
    for i in range(sm.n-1, -1, -1):
        sum = 0
        for idx in range(sm.ai[i],sm.ai[i+1]):
            y[sm.aj[idx]] -= sm.au[idx] * y[i]

def main():
    sm = SparseMatrix()
    inpfile = open("input.txt","r")
    sm.input(inpfile)
    b = list(map(float, inpfile.readline().split()))
    inpfile.close()

    if (not sm.factorizeLU()):
        print('SLAU CANNOT BE SOLVED WITH THIS METHOD')
        exit(0)

    before = time.time()

    find_y(sm, b)
    find_x(sm,b)

    after = time.time()

    print((after-before)*1000.)

    outfile = open("output_py_lin.txt","w")
    for bel in b:
        outfile.write('%.6f\n' % bel)
    outfile.close()

    exit(0)

if __name__ == "__main__":
    main()