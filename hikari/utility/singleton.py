class Singleton(type):
    """
    A metaclass based on https://stackoverflow.com/q/6760685/.
    Declare a singleton class using `MyClass(BaseClass, metaclass=Singleton)`
    """
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
