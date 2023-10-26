class LinkedAttribute:
    """ Class that allows to make the simulation a "reactive" programming software, where the attributes of classes change dynamically.
     One change in a state automatically propagates to other places in the code.
     This is a specialised descriptor class - it will "know" the source of an attribute, and always go to its source to get or set its value.
     In this way the attributes "taken" from objects passed from another class are dynamically updated if the instance of the passed object is modified."""

    def __init__(self, container_path):
        self.container_path = container_path

    def __set_name__(self, owner, name):
        self.name = name

    def _get_container(self, instance):
        container = instance
        for component in self.container_path.split("."):
            container = getattr(container, component)
        return container

    def __get__(self, instance, owner):
        if instance is None:
            return self
        container = self._get_container(instance)
        return getattr(container, self.name)

    def __set__(self, instance, value):
        container = self._get_container(instance)
        setattr(container, self.name, value)

