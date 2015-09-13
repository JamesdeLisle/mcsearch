class a:

    def __init__(self):

        self.a = {'x':1}

trial = a()
print(trial.a['x'])
trial.a['x'] = 0
print(trial.a['x'])
