#!/usr/bin/env python3

class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)
