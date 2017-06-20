# specs.py
"""Volume IB: Testing.
<Name>
<Date>
"""
import math
import sys


# Problem 1 Write unit tests for addition().
# Be sure to install pytest-cov in order to see your code coverage change.
def addition(a,b):
    return a + b

def smallest_factor(n):
    """Finds the smallest prime factor of a number.
    Assume n is a positive integer.
    """
    if n == 1:
        return 1
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return i
    return n


# Problem 2 Write unit tests for operator().
def operator(a, b, oper):
    if type(oper) != str:
        raise ValueError("Oper should be a string")
    if len(oper) != 1:
        raise ValueError("Oper should be one character")
    if oper == "+":
        return a+b
    if oper == "/":
        if b == 0:
            raise ValueError("You can't divide by zero!")
        return a/float(b)
    if oper == "-":
        return a-b
    if oper == "*":
        return a*b
    else:
        raise ValueError("Oper can only be: '+', '/', '-', or '*'")

# Problem 3 Write unit test for this class.
class ComplexNumber(object):
    def __init__(self, real=0, imag=0):
        self.real = real
        self.imag = imag

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def norm(self):
        return math.sqrt(self.real**2 + self.imag**2)

    def __add__(self, other):
        real = self.real + other.real
        imag = self.imag + other.imag
        return ComplexNumber(real, imag)

    def __sub__(self, other):
        real = self.real - other.real
        imag = self.imag - other.imag
        return ComplexNumber(real, imag)

    def __mul__(self, other):
        real = self.real*other.real - self.imag*other.imag
        imag = self.imag*other.real + other.imag*self.real
        return ComplexNumber(real, imag)

    def __div__(self, other):
        if other.real == 0 and other.imag == 0:
            raise ValueError("Cannot divide by zero")
        bottom = (other.conjugate()*other*1.).real
        top = self*other.conjugate()
        return ComplexNumber(top.real / bottom, top.imag / bottom)

    def __eq__(self, other):
        return self.imag == other.imag and self.real == other.real

    def __str__(self):
        return "{}{}{}i".format(self.real, '+' if self.imag >= 0 else '-',
                                                                abs(self.imag))

# Problem 5: Write code for the Set game here
def set_game(filename):
    possible_vals=['0','1','2']
    with open(filename, 'r') as infile:
        lines = infile.readlines()
        if len(lines) != 12:
            raise ValueError("There must only be 12 cards!")
    duplicates=set()
    for i in xrange(12):
        duplicates.add(lines[i].strip())
    if len(duplicates) != 12:
        raise ValueError("You have duplicate cards.")
    cards = []
    for line in lines:
        cards.append(list(line.strip()))
    for card in cards:
        if len(card) != 4:
            raise ValueError("Each card can only have 4 attributes!")
        else:
            for i in xrange(4):
                if card[i] not in possible_vals:
                    raise ValueError("There can only be 0, 1, or 2!")
                card[i] = int(card[i])
    matches = 0
    for i in xrange(10):
        for j in xrange(i+1,11):
            for k in xrange(j+1,12):
                match = True
                if (cards[i][0]+cards[j][0]+cards[k][0]) %3 != 0:
                    match = False
                    print "1"
                if match and (cards[i][1]+cards[j][1]+cards[k][1])%3!=0:
                    match = False
                    print "2"
                if match and (cards[i][2]+cards[j][2]+cards[k][2])%3!=0:
                    match = False
                    print "3"
                if match and (cards[i][3]+cards[j][3]+cards[k][3])%3!=0:
                    match = False  
                    print "4"
                if match:
                    matches += 1
    return matches
        