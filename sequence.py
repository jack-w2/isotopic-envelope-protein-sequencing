class Seq:
    def __init__(self, seq, type):
        self.seq = seq
        self.type = type

    def __str__(self):
        return f'Seq({self.seq}, {self.type})'

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        if isinstance(other, str):
            return Seq(self.seq+other, self.type)
        assert self.type == other.type
        return Seq(self.seq + other.seq, type = self.type)
