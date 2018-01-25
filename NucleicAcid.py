class IncorrectSequence(Exception):
    """
    Exception raised when input contain forbidden signs
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class SequenceTypeError(Exception):
    """
    Exception raised when nucleic acid type can not be defined
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class NucleicAcid(object):
    """
    Class reprezenting nucleic acid sequence
    """

    DNA = ['A', 'T', 'C', 'G']
    RNA = ['A', 'T', 'C', 'G' 'U']  # T will be changed to U
    MIXED = []  # to be added

    def __init__(self, seq, seq_type):
        self.seq = self.seq_check(seq)
        self.type = self.type_check(seq_type)
        self.length = len(self.seq)

    def seq_check(self, seq):
        if type(seq) == str:
            seq = seq.replace(' ', '')
        else:
            raise TypeError("Sequence must be of type string")

        if seq.isalpha():
            return seq.upper()
        else:
            raise IncorrectSequence("Sequence contain incorrect signs")

    def type_check(self, seq_type):
        if type(seq_type) == str:
            seq_type = seq_type.replace(' ', '').upper()
        else:
            raise TypeError("Sequence must be of type string")

        if self.seq.strip(''.join(self.DNA)) == '' and seq_type == "DNA":
            return seq_type
        elif self.seq.strip(''.join(self.RNA)) == '' and seq_type == "RNA":
            self.seq = self.seq.replace('T', 'U')
            return seq_type
        else:
            raise SequenceTypeError("Type of the sequence cannot be confirmed \
                                    only DNA and RNA are allowed.")


if __name__ == "__main__":
    pass
