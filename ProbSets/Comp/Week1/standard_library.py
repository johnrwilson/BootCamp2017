# standard_library.py
"""Introductory Labs: The Standard Library
<John>
<Math 345>
<9/1/16>
"""

import box

from sys import argv

from random import randint

import calculator as calc

# Problem 1
def prob1(l):
    """Accept a list 'l' of numbers as input and return a new list with the
    minimum, maximum, and average of the contents of 'l'.
    """
    output = [min(l), max(l), float(sum(l))/len(l)]
    return output


# Problem 2
def prob2():
    """Programmatically determine which Python objects are mutable and which
    are immutable. Test numbers, strings, lists, tuples, and dictionaries.
    Print your results to the terminal.
    """
    num1 = 8
    num2 = num1
    num2 += 1
    print "Numbers are mutable: ", num1 == num2
    
    word1 = "Car"
    word2 = word1
    word2 += 'a'
    print "Strings are mutable: ", word1 == word2

    list1 = ["a", "b", "c", "d", "e", "f"]
    list2 = list1
    list2.append(1)
    print "Lists are mutable: ", list1 == list2

    tup1 = ('green', 'yellow', 'orange')
    tup2 = tup1
    tup2 += (1,)
    print "Tuples are mutable: ", tup1 == tup2

    dict1 = {1: 'x', 2: 'b'}
    dict2 = dict1
    dict2[1] = 'a'
    print "Dictionaries are mutable: ", dict1 == dict2


# Problem 3: Create a 'calculator' module and implement this function.
def prob3(a,b):
    """Calculate and return the length of the hypotenuse of a right triangle.
    Do not use any methods other than those that are imported from your
    'calculator' module.

    Parameters:
        a (float): the length one of the sides of the triangle.
        b (float): the length the other nonhypotenuse side of the triangle.

    Returns:
        The length of the triangle's hypotenuse.
    """
    return calc.sqrt(calc.add(calc.multiply(a, a), calc.multiply(b, b)))


# Problem 4: Implement shut the box.

if __name__ == "__main__":

    if len(argv) == 2:
        script, user_name = argv
    else:
        script = argv
        user_name = raw_input("What is your name? ")

    is_valid = True
    numbers = range(1, 10)

    while is_valid:
        print "\nNumbers left: ", numbers
        remaining = sum(numbers)
        if remaining > 6:
            roll = randint(1, 6) + randint(1, 6)
        else:
            roll = randint(1, 6)
        print "Roll: ", roll
        is_valid = box.isvalid(roll, numbers)
        if not is_valid:
            print "Score for player", user_name, ":", sum(numbers)
            if sum(numbers) == 0:
                print "Congratulations, you shut the box!"		
        else:
            nums_to_delete = box.parse_input(raw_input("Numbers to eliminate: "), numbers)
            while sum(nums_to_delete) != roll:
                print "Invalid input, try again!"
                nums_to_delete = box.parse_input(raw_input("Numbers to eliminate: "), numbers)
            else:
                for n in nums_to_delete:
                    numbers.remove(n)