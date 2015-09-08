import random

class coeffContainer:

    def __init__(self,real,imag,flag_real,flag_imag):

        self.imag = imag
        self.real = real
        self.imag_flag = flag_imag
        self.real_flag = flag_real
        self.value = self.real + 1j * self.imag

def getUniformSample(flag,width):

    if flag:
        return random.uniform(-width,width)    
    else:
        return 1
   
def getGaussianSample(flag,standard_deviation,current_value):

    if flag:
        return random.gauss(current_value,standard_deviation)
    else:
        return 1

class Lambda:

    def __init__(self,coefficient_names):

        self.coefficients = {}
        for coeff in coefficient_names:
            flag_real = coefficient_names[coeff][0]
            flag_imag = coefficient_names[coeff][1]
            self.coefficients[coeff] = coeffContainer(getUniformSample(flag_real,10.0),getUniformSample(flag_imag,10.0),flag_real,flag_imag)

    def walk(self,standard_deviation):

        for coeff in self.coefficients:
            self.coefficients[coeff].real = getGaussianSample(self.coefficients[coeff].real_flag,standard_deviation,self.coefficients[coeff].real)
            self.coefficients[coeff].imag = getGaussianSample(self.coefficients[coeff].imag_flag,standard_deviation,self.coefficients[coeff].imag)


        


