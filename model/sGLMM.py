import scipy.optimize as opt
import time

from ProximalGradientDescent import ProximalGradientDescent
from GroupLasso import GFlasso

from helpingMethods import *
def normalize(x):
    m = np.mean(x)
    s = np.std(x)
    return (x - m) / s

class sGLMM:
    def __init__(self, numintervals=100, ldeltamin=-5, ldeltamax=5, discoverNum=None, scale=0, learningRate=1e-5,
                 lam=1, reg_min=1e-7, reg_max=1e7, threshold=0.618, isQuiet=False, cv_flag=False):
        self.numintervals = numintervals
        self.ldeltamin = ldeltamin
        self.ldeltamax = ldeltamax
        self.discoverNum = discoverNum
        self.scale = scale

        self.learningRate = learningRate
        self.lam = lam
        self.reg_min = reg_min
        self.reg_max = reg_max
        self.threshold = threshold
        self.isQuiet = isQuiet
        self.cv_flag = cv_flag
        self.beta = None

    def cross_val_score(self, X, y, lam, cv=5):
        scores = []
        [n, p] = X.shape
        b = n / cv
        for i in range(cv):
            if not self.isQuiet: print "\tRun", i, "..."
            ind = np.arange(b) + b * i
            Xtr = np.delete(X, ind, axis=0)
            ytr = np.delete(y, ind, axis=0)
            Xte = X[ind, :]
            yte = y[ind]
            self.beta = self.runLasso(X=Xtr, Y=ytr, lam_=lam)
            ypr = self.predict(Xte)
            s = np.mean(np.square(ypr - yte))
            scores.append(s)
        return scores

    def crossValidation(self, X, y):
        minError = np.inf
        minLam = 0
        for i in range(-8, 8):
            lam = np.power(10., i)
            # model.set_lambda(lam)
            if not self.isQuiet: print "lambda for this iteration: ", lam
            scores = self.cross_val_score(X, y, lam, cv=5)
            score = np.mean(np.abs(scores))
            if not self.isQuiet: print "score: ", score
            if score < minError:
                minError = score
                minLam = lam
        print "The best lambda is ", minLam
        beta = self.runLasso(X=X, Y=y, lam_=minLam)
        return beta

    def train(self, X, K, y):
        print 'Computing ...'
        Kva, Kve = np.linalg.eigh(K)
        time_start = time.time()
        [n_s, n_f] = X.shape
        assert X.shape[0] == y.shape[0], 'dimensions do not match'
        assert K.shape[0] == K.shape[1], 'dimensions do not match'
        assert K.shape[0] == X.shape[0], 'dimensions do not match'
        if y.ndim == 1:
            y = scipy.reshape(y, (n_s, 1))

        S, U, ldelta0, monitor_nm = self.train_nullmodel(y, K, S=Kva, U=Kve)

        delta0 = scipy.exp(ldelta0)
        Sdi = 1. / (S + delta0)
        Sdi_sqrt = scipy.sqrt(Sdi)
        SUX = scipy.dot(U.T, X)
        SUX = SUX * scipy.tile(Sdi_sqrt, (n_f, 1)).T
        SUy = scipy.dot(U.T, y)
        SUy = SUy * scipy.reshape(Sdi_sqrt, (n_s, 1))
        # the normalize part can be commented
        for i in range(0, y.shape[1]):
            SUy[:, i]=normalize(SUy[:, i])
        SUX = normalize(SUX)

        if self.cv_flag:
            self.beta=self.crossValidation(SUX, SUy)
        elif self.discoverNum is not None:
            self.beta = self.cv_train(X=SUX, Y=SUy, K=int(self.discoverNum)*y.shape[1])
        else:
            self.beta = self.runLasso(SUX, SUy, lam_=self.lam)

        time_end = time.time()
        time_diff = time_end - time_start

        print 'Computation completed in %.2fs' % (time_diff)
        return self.beta

    def cv_train(self, X, Y, K):
        regMin = self.reg_min
        regMax = self.reg_max
        betaM = None
        iteration = 0
        patience = 100
        ss = []
        time_start = time.time()
        time_diffs = []
        minFactor = 0.5
        maxFactor = 2
        while regMin + 1e-5 < regMax and iteration < patience:
            iteration += 1
            reg = np.exp((np.log(regMin)+np.log(regMax)) / 2.0)
            coef_=self.runLasso(X,Y, reg)
            k = len(np.where(coef_ != 0)[0])
            if not self.isQuiet: print "\tIter:%d\t   lambda:%.5f  non-zeroes:%d" % (iteration, reg, k)
            ss.append((reg, k))
            if k < K * minFactor:    # Regularizer too strong
                if not self.isQuiet: print '\tRegularizer is too strong, shrink the lambda'
                regMax = reg
            elif k > K * maxFactor:  # Regularizer too weak
                if not self.isQuiet: print '\tRegularization is too weak, enlarge lambda'
                regMin = reg
            else:
                betaM = coef_
                break
        time_diffs.append(time.time() - time_start)
        time_start = time.time()
        return betaM


    def runLasso(self, X, Y,lam_):
        pgd = ProximalGradientDescent(learningRate=self.learningRate * 2e4)
        model = GFlasso(lambda_flasso=lam_, gamma_flasso=0.7, mau=0.1)

        # Set X, Y, correlation
        model.setXY(X, Y)
        graph_temp = np.cov(Y.T)
        graph_temp = graph_temp.reshape((Y.shape[1], Y.shape[1]))
        graph = np.zeros((Y.shape[1], Y.shape[1]))
        for i in range(0, Y.shape[1]):
            for j in range(0, Y.shape[1]):
                graph[i, j] = graph_temp[i, j] / (np.sqrt(graph_temp[i, i]) * (np.sqrt(graph_temp[j, j])))
                if (graph[i, j] < self.threshold):
                    graph[i, j] = 0
        model.corr_coff = graph

        pgd.run(model)
        return model.beta


    def hypothesisTest(self, UX, Uy, X, UX0, X0):
        [m, n] = X.shape
        p = []
        for i in range(n):
            if UX0 is not None:
                UXi = np.hstack([UX0 ,UX[:, i].reshape(m, 1)])
                XX = matrixMult(UXi.T, UXi)
                XX_i = linalg.pinv(XX)
                beta = matrixMult(matrixMult(XX_i, UXi.T), Uy)
                Uyr = Uy - matrixMult(UXi, beta)
                Q = np.dot( Uyr.T, Uyr)
                sigma = Q * 1.0 / m
            else:
                Xi = np.hstack([X0 ,UX[:, i].reshape(m, 1)])
                XX = matrixMult(Xi.T, Xi)
                XX_i = linalg.pinv(XX)
                beta = matrixMult(matrixMult(XX_i, Xi.T), Uy)
                Uyr = Uy - matrixMult(Xi, beta)
                Q = np.dot(Uyr.T, Uyr)
                sigma = Q * 1.0 / m
            ts, ps = tstat(beta[1], XX_i[1, 1], sigma, 1, m)
            if -1e10 < ts < 1e10:
                p.append(ps)
            else:
                p.append(1)
        return p

    def predict(self, X):
        return np.dot(X, self.beta)

    def train_nullmodel(self, y, K, S=None, U=None):
        self.ldeltamin += self.scale
        self.ldeltamax += self.scale

        if S is None or U is None:
            S, U = linalg.eigh(K)

        Uy = scipy.dot(U.T, y)
        # grid search
        nllgrid = scipy.ones(self.numintervals + 1) * scipy.inf
        ldeltagrid = scipy.arange(self.numintervals + 1) / (self.numintervals * 1.0) * (self.ldeltamax - self.ldeltamin) + self.ldeltamin
        for i in scipy.arange(self.numintervals + 1):
            nllgrid[i] = nLLeval(ldeltagrid[i], Uy, S)

        nllmin = nllgrid.min()
        ldeltaopt_glob = ldeltagrid[nllgrid.argmin()]

        for i in scipy.arange(self.numintervals - 1) + 1:
            if (nllgrid[i] < nllgrid[i - 1] and nllgrid[i] < nllgrid[i + 1]):
                ldeltaopt, nllopt, iter, funcalls = opt.brent(nLLeval, (Uy, S),
                                                              (ldeltagrid[i - 1], ldeltagrid[i], ldeltagrid[i + 1]),
                                                              full_output=True)
                if nllopt < nllmin:
                    nllmin = nllopt
                    ldeltaopt_glob = ldeltaopt


        return S, U, ldeltaopt_glob, None
