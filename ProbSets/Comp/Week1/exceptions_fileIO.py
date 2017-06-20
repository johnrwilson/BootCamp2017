# exceptions_fileIO.py
"""Introductory Labs: Exceptions and File I/O.
<John Wilson>
<Math 345>
<9/5/2016>
"""

from random import choice

# Problem 1
def arithmagic():
    step_1 = raw_input("Enter a 3-digit number where the first and last "
                                            "digits differ by 2 or more: ")
    if (len(step_1) != 3):
        raise ValueError("The number must be 3 digits")
    if abs(int(step_1[0])-int(step_1[2])) < 2:
        raise ValueError("The first and last digits must differ by 2 or more")
        
    step_2 = raw_input("Enter the reverse of the first number, obtained "
                                            "by reading it backwards: ")
    if step_2 != step_1[::-1]:
        raise ValueError("The number entered must be the reverse of %s" % (step_1))
        
    step_3 = raw_input("Enter the positive difference of these numbers: ")
    if int(step_3)!= abs(int(step_2)-int(step_1)):
        raise ValueError("The number must be the positive difference between the first two number")
        
    step_4 = raw_input("Enter the reverse of the previous result: ")
    if step_4 != step_3[::-1]:
        raise ValueError("The number entered must be the reverse of %s" % (step_3))
      
    print str(step_3) + " + " + str(step_4) + " = 1089 (ta-da!)"


# Problem 2
def random_walk(max_iters=1e9):
    walk = 0
    direction = [-1, 1]
    try:
        for i in xrange(int(max_iters)):
            walk += choice(direction)
    except KeyboardInterrupt:
        print "Process interrupted at iteration %d" % (i)
    else:
        print "Process completed."
    finally:
        return walk

        
# Problems 3 and 4: Write a 'ContentFilter' class.
class ContentFilter(object):
    def __init__(self, file_name):
        if not isinstance(file_name, str):
            raise TypeError("The argument must be a string!")
        self.file = file_name
        try:
            myfile = open(self.file, 'r')
            self.contents = myfile.read()
        finally:
            myfile.close()
            
    def uniform(self, name, mode='w', case="upper"):
        if mode != 'w' and mode != 'a':
            raise ValueError("The mode entered must be either \'w\' or \'a\'.")
        self.case = case
        if self.case != "upper" and self.case!= "lower":
            raise ValueError("You must enter either \"upper\" or \"lower\" as your third argument")
        else:
            if self.case == "lower":
                text = self.contents.lower()
            elif self.case == "upper":
                text = self.contents.upper()
            with open(name, mode) as outfile:
                outfile.write(text)
    
    def reverse(self, name, mode='w', unit="line"):
        if mode != 'w' and mode != 'a':
            raise ValueError("The mode entered must be either \'w\' or \'a\'.")
        if unit != "line" and unit!= "word":
            raise ValueError("You must enter either \"line\" or \"word\" as your third argument")
        else:
            if unit == "line":
                text = self.contents.split('\n')
                text.reverse()
                result = "\n".join(text)
            elif unit == "word":
                text = self.contents.split('\n')
                for i in xrange(len(text)):
                    temp = text[i].split()
                    temp.reverse()
                    text[i] = ' '.join(temp)
                result = "\n".join(text)
            with open(name, mode) as outfile:
                outfile.write(result)
        
    def transpose(self, name, mode='w'):
        if mode != 'w' and mode != 'a':
            raise ValueError("The mode entered must be either \'w\' or \'a\'.")
        else:
            lines = self.contents.split('\n')
            max_words=0
            final_lines=[]
            for i in xrange(len(lines)):
                temp = lines[i].split()
                if len(temp) > max_words:
                    max_words = len(temp)
            for j in xrange(max_words):
                words=[]
                for i in xrange(len(lines)):
                    temp = lines[i].split()
                    if j<len(temp):
                        words.append(temp[j])
                final_lines.append(' '.join(words))
            result = '\n'.join(final_lines)
            with open(name, mode) as outfile:
                outfile.write(result)
                
    def __str__(self):
        characters = len(self.contents)
        
        letter=[]
        for i in xrange(characters):
            if self.contents[i].isalpha():
                letter.append(self.contents[i])
        letters = len(letter)
        
        number=[]
        for i in xrange(characters):
            if self.contents[i].isdigit():
                number.append(self.contents[i])
        nums = len(number)
        
        whitespace = self.contents.count('\n') + self.contents.count('\t') + self.contents.count(' ')
        
        lines = len(self.contents.split('\n'))
        
        return "Source file: \t\t%s\nTotal characters: \t%d\nAlphabetic characters: \t%d\nNumerical characters: \t%d\nWhitespace characters: \t%d\nNumber of lines: \t%d" % (self.file,characters,letters,nums,whitespace,lines)
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        