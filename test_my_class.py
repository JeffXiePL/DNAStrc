from my_class import Person
import pytest

# To run these tests, from the command line, run:
# pytest test_my_class.py


def test_assignment():
    person = Person("Alice", 30)
    assert person.name == "Alice"
    assert person.age == 30


def test_str():
    person = Person("Alice", 30)
    assert str(person) == "Alice is 30 years old."


def test_validate_types():
    with pytest.raises(TypeError):
        Person(1, 30)
    with pytest.raises(TypeError):
        Person("Alice", "30")
