class Person:
    def __init__(self, name, surname) -> None:
        self.name = name
        self.surname = surname

    def from_string(input: str):
        name, surname = input.split(" ")
        return print(Person(name, surname))
    
    @classmethod
    def from_string2(cls, input):
        name, surname = input.split(" ")
        return cls(name, surname)

person = Person.from_string("Peilun Xie")
