# -*- coding: utf-8 -*-
"""
Simple progress bar that uses only Python Standard Library.
Code obtained from Gnu/Linux Magazine HS 86 Mémo Python.
"""
import time

class ProgressBar:
    def __init__(self, msg = "Work in progress"):
        self.current = 0
        self.msg = msg
    
    def add(self, increment):
        self.current += increment
        if self.current > 100:
            self.current = 100
    
    def show(self):
        text = "{} : [{}%]".format(self.msg, self.current)
        if self.isFinished():
            text += "\n"
        else:
            text += "\r"
        print(text, end="")
        
    def isFinished(self):
        return self.current == 100
    
if __name__ == "__main__":
    bar = ProgressBar()
    
    while not bar.isFinished():
        bar.show()
        time.sleep(1)
        bar.add(5)