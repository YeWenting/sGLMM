__author__ = "Xiang Liu"


class ProximalGradientDescent:
    def __init__(self, tolerance=0.000001, learningRate=0.001, prev_residue=999999999999.,
                 maxIteration=1000):
        self.tolerance = tolerance
        self.learningRate = learningRate
        self.prev_residue = prev_residue
        self.maxIteration = maxIteration

    def stop(self):
        self.shouldStop = True

    def setUpRun(self):
        self.isRunning = True
        self.progress = 0.
        self.shouldStop = False

    def finishRun(self):
        self.isRunning = False
        self.progress = 1.

    def run(self, model):
        epoch = 0
        self.learningRate = self.learningRate *1 #5e3
        residue = model.cost()
        theta = 1.
        theta_new = 0.
        beta_prev = model.beta
        beta_curr = model.beta
        beta = model.beta
        beta_best = model.beta
        diff = self.tolerance * 2
        residue_best = 9999999999
        x = []
        y = []
        while (epoch < self.maxIteration and diff > self.tolerance):
            epoch = epoch + 1
            theta_new = 2. / (epoch + 3.)
            grad = model.gradient()
            in_ = beta - 1. / model.getL() * grad
            beta_curr = model.proximal_operator(in_, self.learningRate)
            beta = beta_curr + (1 - theta) / theta * theta_new * (beta_curr - beta_prev)
            beta_prev = beta_curr
            theta = theta_new
            model.beta = beta
            residue = model.cost()
            diff = abs(self.prev_residue - residue)
            self.prev_residue = residue
            # print epoch, residue
            # print beta
            if (residue < residue_best):
                beta_best = beta
                residue_best = residue
            # roc(model.beta, str)
            x.append(epoch)
            y.append(residue)
        model.beta = beta_best
        return residue_best