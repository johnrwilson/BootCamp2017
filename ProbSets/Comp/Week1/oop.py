# oop.py
"""Introductory Labs: Object Oriented Programming.
<John Wilson>
<Math 321 Section 1>
<8/30/16>
"""

class Backpack(object):
    """A Backpack object class. Has a name and a list of contents.

    Attributes:
        name (str): the name of the backpack's owner.
        contents (list): the contents of the backpack.
    """

    # Problem 1: Modify __init__() and put(), and write dump().
    def __init__(self, name, color, max_size=5):
        """Set the name and initialize an empty contents list.

        Inputs:
            name (str): the name of the backpack's owner.

        Returns:
            A Backpack object wth no contents.
        """
        self.name = name
        self.contents = []
        self.color = color
        self.max_size = max_size

    def put(self, item):
        """Add 'item' to the backpack's list of contents."""
        if len(self.contents) < self.max_size:
            self.contents.append(item)
        else:
            print "No Room!"

    def take(self, item):
        """Remove 'item' from the backpack's list of contents."""
        self.contents.remove(item)
        
    def dump(self):
        self.contents = []
    
        

    # Magic Methods -----------------------------------------------------------

    # Problem 3: Write __eq__() and __str__().
    def __add__(self, other):
        """Add the number of contents of each Backpack."""
        return len(self.contents) + len(other.contents)

    def __lt__(self, other):
        """Compare two backpacks. If 'self' has fewer contents
        than 'other', return True. Otherwise, return False.
        """
        return len(self.contents) < len(other.contents)
        
    def __eq__(self, other):
        """Compare two backpacks. Returns True if both backpacks
        have the same number of contents.
        """
        return len(self.contents) == len(other.contents) and self.name == other.name and self.color == other.color
        
    def __str__(self):
        return "Owner:\t\t%s\nColor:\t\t%s\nSize:\t\t%d\nContents:\t[<" % (self.name, self.color, self.max_size) + ">, <".join(self.contents) + ">]"
        


# An example of inheritance. You are not required to modify this class.
class Knapsack(Backpack):
    """A Knapsack object class. Inherits from the Backpack class.
    A knapsack is smaller than a backpack and can be tied closed.

    Attributes:
        name (str): the name of the knapsack's owner.
        color (str): the color of the knapsack.
        max_size (int): the maximum number of items that can fit
            in the knapsack.
        contents (list): the contents of the backpack.
        closed (bool): whether or not the knapsack is tied shut.
    """
    def __init__(self, name, color, max_size=3):
        """Use the Backpack constructor to initialize the name, color,
        and max_size attributes. A knapsack only holds 3 item by default
        instead of 5.

        Inputs:
            name (str): the name of the knapsack's owner.
            color (str): the color of the knapsack.
            max_size (int): the maximum number of items that can fit
                in the knapsack. Defaults to 3.

        Returns:
            A Knapsack object with no contents.
        """
        Backpack.__init__(self, name, color, max_size)
        self.closed = True

    def put(self, item):
        """If the knapsack is untied, use the Backpack.put() method."""
        if self.closed:
            print("I'm closed!")
        else:
            Backpack.put(self, item)

    def take(self, item):
        """If the knapsack is untied, use the Backpack.take() method."""
        if self.closed:
            print("I'm closed!")
        else:
            Backpack.take(self, item)


# Problem 2: Write a 'Jetpack' class that inherits from the 'Backpack' class.
class Jetpack(Backpack):
    def __init__(self, name, color, max_size=2, fuel = 10):
        """Initializes the jetpack as a backpack with a fuel attribute
        """
        Backpack.__init__(self, name, color, max_size)
        self.fuel = fuel
        
    def fly(self, fuel):
        """Subtracts the amount of fuel used, or if there is less fuel than
        passed to the function, returns error message
        """
        if fuel < self.fuel:
            self.fuel -= fuel
        else:
            print "Not enough fuel!"
    
    def dump(self):
        Backpack.dump(self)
        self.fuel = 0

# Problem 4: Write a 'ComplexNumber' class.
class ComplexNumber(object):
    """A class that stores both the real and imaginary parts of a complex number
    """
    def __init__(self,real=0,imag=0):
        """Initializes the complex number, with the first argument saved as the 
        real component, and the second argument saved as the imaginary component
        """
        self.real = real
        self.imag = imag
        
    def conjugate(self):
        """Returns the complex conjugate of the complex number
        """
        return ComplexNumber(self.real, -1 * self.imag)
        
    def __abs__(self):
        """Returns the magnitude of the complex number, obtained by taking the 
        square root of the real component squared added to the imaginary component
        squared
        """
        return (self.real**2 + self.imag**2)**.5
        
    def __lt__(self, other):
        """Returns true if the complex conjugate of the first complex number
        is less than the complex conjugate of the second complex number
        """
        return (abs(self) < abs(other))
        
    def __gt__(self, other):
        """Returns true if the complex conjugate of the first complex number
        is greater than the complex conjugate of the second complex number
        """
        return (abs(self) > abs(other))
        
    def __eq__(self, other):
        """Returns true if the complex numbers entered are equal
        """
        return (self.real == other.real and self.imag == other.imag)
        
    def __ne__(self, other):
        """Returns true if the complex numbers entered are not equal
        """
        return not ComplexNumber.__eq__(self, other)
    
    def __add__(self, other):
        """Returns the sum of the two complex numbers
        """
        return ComplexNumber(self.real+other.real, self.imag+other.imag)
        
    def __sub__(self, other):
        """Returns the difference of the compley numbers
        """
        return ComplexNumber(self.real-other.real, self.imag-other.imag)
        
    def __mul__(self, other):
        """Returns the product of the two imaginary numbers
        """
        real_part = (self.real*other.real) - (self.imag*other.imag)
        imag_part = (self.real*other.imag) + (self.imag*other.real)
        return ComplexNumber(real_part, imag_part)
        
    def __div__(self, other):
        """Returns the quotient of the two imaginary numbers
        """
        den = other.real**2 + other.imag**2
        real_part = ((self.real*other.real)+(self.imag*other.imag))/float(den)
        imag_part = (self.imag*other.real-self.real*other.imag)/float(den)
        return ComplexNumber(real_part, imag_part)

    #the built-in Python complex number class has the advantage of resistering as a complex number if you just put a j after an int.
    