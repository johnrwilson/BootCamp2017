# test_specs.py
"""Volume 1B: Testing.
<John Wilson>
<Math 347>
<1/10/17>
"""

import specs
import pytest
import os

# Problem 1: Test the addition and fibonacci functions from specs.py
def test_addition():
    assert 1 + 2 == specs.addition(1,2)
    assert 4 + -10 == specs.addition(4,-10)

def test_smallest_factor():
    assert specs.smallest_factor(8) == 2
    assert specs.smallest_factor(1) == 1
    assert specs.smallest_factor(7) == 7

# Problem 2: Test the operator function from specs.py
def test_operator():
    assert specs.operator(1,2,'+') == 3
    assert specs.operator(3,1,'/') == 3
    assert specs.operator(3,1,'*') == 3
    assert specs.operator(3,1,'-') == 2
    with pytest.raises(Exception) as excinfo1:
        specs.operator(2,3,4)
    with pytest.raises(Exception) as excinfo2:
        specs.operator(4,3,'+-')
    with pytest.raises(Exception) as excinfo3:
        specs.operator(3,0,'/')
    with pytest.raises(Exception) as excinfo4:
        specs.operator(4,2,'a')
    assert excinfo1.typename == 'ValueError'
    assert excinfo2.typename == 'ValueError'
    assert excinfo3.typename == 'ValueError'
    assert excinfo4.typename == 'ValueError'
    assert excinfo1.value.args[0] == "Oper should be a string"
    assert excinfo2.value.args[0] == "Oper should be one character"
    assert excinfo3.value.args[0] == "You can't divide by zero!"
    assert excinfo4.value.args[0] == "Oper can only be: '+', '/', '-', or '*'"

# Problem 3: Finish testing the complex number class
@pytest.fixture
def set_up_complex_nums():
    number_1 = specs.ComplexNumber(1, 2)
    number_2 = specs.ComplexNumber(5, 5)
    number_3 = specs.ComplexNumber(2, 9)
    return number_1, number_2, number_3

def test_complex_addition(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1 + number_2 == specs.ComplexNumber(6, 7)
    assert number_1 + number_3 == specs.ComplexNumber(3, 11)
    assert number_2 + number_3 == specs.ComplexNumber(7, 14)
    assert number_3 + number_3 == specs.ComplexNumber(4, 18)

def test_complex_multiplication(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1 * number_2 == specs.ComplexNumber(-5, 15)
    assert number_1 * number_3 == specs.ComplexNumber(-16, 13)
    assert number_2 * number_3 == specs.ComplexNumber(-35, 55)
    assert number_3 * number_3 == specs.ComplexNumber(-77, 36)

def test_complex_subtraction(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1 - number_2 == specs.ComplexNumber(-4, -3)

def test_complex_conjugate(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1.conjugate() == specs.ComplexNumber(1, -2)

def test_complex_norm(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1.norm() == (number_1.real**2 + number_1.imag**2)**.5

def test_complex_equals(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    num = specs.ComplexNumber(1, 2)
    assert number_1 == num

def test_complex_to_string(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    assert number_1.__str__() == '1+2i'

def test_complex_division(set_up_complex_nums):
    number_1, number_2, number_3 = set_up_complex_nums
    num = number_1
    num2 = specs.ComplexNumber(0, 0)
    with pytest.raises(Exception) as excinfo:
        number_1.__div__(num2)
    assert number_1/num == specs.ComplexNumber(1.0,0)  
    assert excinfo.typename == 'ValueError'
    assert excinfo.value.args[0] == "Cannot divide by zero"
    
# Problem 4: Write test cases for the Set game.
def test_set_name():
    
    with pytest.raises(Exception) as excinfo:
        specs.set_game('notarealfile.txt')
    assert excinfo.typename == 'IOError'
    #assert excinfo.value.args[0] == "[Errno 2] No such file or directory: 'notarealfile.txt'"
    
    with pytest.raises(Exception) as excinfo:
        specs.set_game('hands/duplicates.txt')
    assert excinfo.typename == 'ValueError'
    assert excinfo.value.args[0] == "You have duplicate cards."

    with pytest.raises(Exception) as excinfo:
        specs.set_game('hands/11_cards.txt')
    assert excinfo.typename == 'ValueError'
    assert excinfo.value.args[0] == "There must only be 12 cards!"
    
    with pytest.raises(Exception) as excinfo:
        specs.set_game('hands/invalid_numbers.txt')
    assert excinfo.typename == 'ValueError'
    assert excinfo.value.args[0] == "There can only be 0, 1, or 2!"
    
    with pytest.raises(Exception) as excinfo:
        specs.set_game('hands/bad_cards.txt')
    assert excinfo.typename == 'ValueError'
    assert excinfo.value.args[0] == "Each card can only have 4 attributes!"
    
    assert specs.set_game('hands/set.txt') == 6