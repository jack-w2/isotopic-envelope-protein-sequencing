class Seq:
    def __init__(self, seq, type):
        self.seq = seq
        self.type = type

    def __str__(self):
        return f'{self.seq}, {self.type}'


