class Person:
    def __init__(self, name: str, age: int):
        if not isinstance(name, str):
            raise TypeError("Name must be a string.")
        if not isinstance(age, int):
            raise TypeError("Age must be an integer.")
        self.name = name
        self.age = age

    def __str__(self):
        return f"{self.name} is {self.age} years old."
