import functools

class FunctionalExpression:
    def __init__(self, parts):
        self.parts = tuple(parts)
    def dump(self, indent="  "):
        print("%s%s" % (indent, self._dump()))
        for part in self.parts:
            part.dump(indent + "  ")
    def _dump(self):
        return self.__class__.__name__
    def instantiate(self, var_mapping, init_facts):
        raise ValueError("Cannot instantiate condition: not normalized")

class NumericConstant(FunctionalExpression):
    parts = ()
    def __init__(self, value):
        self.value = float(value)

    def __hash__(self):
        return hash(self.value)
    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self.value == other.value)
    def __str__(self):
        return "%s %s" % (self.__class__.__name__, self.value)
    def _dump(self):
        return str(self)
    def instantiate(self, var_mapping, init_facts):
        return self


class PrimitiveNumericExpression(FunctionalExpression):
    parts = ()
    arithmetic_operations = {
        "*": lambda args: functools.reduce(lambda a, b: a * b, args),
        "/": lambda args: functools.reduce(lambda a, b: a / b, args),
        "+": lambda args: functools.reduce(lambda a, b: a + b, args),
        "-": lambda args: functools.reduce(lambda a, b: a - b, args)
    }

    def __init__(self, symbol, args):
        self.symbol = symbol

        if self.symbol in self.arithmetic_operations.keys():
            self.args = tuple([NumericConstant(float(a)) if not isinstance(a, list) else PrimitiveNumericExpression(a[0], a[1:]) for a in args])
        else:
            self.args = tuple(args)
        self.hash = hash((self.__class__, self.symbol, self.args))
    def __hash__(self):
        return self.hash
    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self.symbol == other.symbol
                and self.args == other.args)
    def __str__(self):
        return "%s %s(%s)" % ("PNE", self.symbol, ", ".join(map(str, self.args)))
    def dump(self, indent="  "):
        print("%s%s" % (indent, self._dump()))
    def _dump(self):
        return str(self)
    def instantiate(self, var_mapping, init_assignments):

        def getvalue(exp):
            if isinstance(exp, NumericConstant):
                return exp.value
            else:
                return exp.instantiate(var_mapping, init_assignments).value

        if self.symbol in self.arithmetic_operations.keys():
            op = self.arithmetic_operations[self.symbol]
            return NumericConstant(op([getvalue(a) for a in self.args]))
        else:
            args = [var_mapping.get(arg, arg) for arg in self.args]
            pne = PrimitiveNumericExpression(self.symbol, args)
            assert self.symbol != "total-cost"
            # We know this expression is constant. Substitute it by corresponding
            # initialization from task.
            result = init_assignments.get(pne)
            assert result is not None, "Could not find instantiation for PNE: %r" % (str(pne),)
            return result

class FunctionAssignment:
    def __init__(self, fluent, expression):
        self.fluent = fluent
        self.expression = expression
    def __str__(self):
        return "%s %s %s" % (self.__class__.__name__, self.fluent, self.expression)
    def dump(self, indent="  "):
        print("%s%s" % (indent, self._dump()))
        self.fluent.dump(indent + "  ")
        self.expression.dump(indent + "  ")
    def _dump(self):
        return self.__class__.__name__
    def instantiate(self, var_mapping, init_facts):
        if not (isinstance(self.expression, PrimitiveNumericExpression) or
                isinstance(self.expression, NumericConstant)):
            raise ValueError("Cannot instantiate assignment: not normalized")
        # We know that this assignment is a cost effect of an action (for initial state
        # assignments, "instantiate" is not called). Hence, we know that the fluent is
        # the 0-ary "total-cost" which does not need to be instantiated
        assert self.fluent.symbol == "total-cost"
        fluent = self.fluent
        expression = self.expression.instantiate(var_mapping, init_facts)
        return self.__class__(fluent, expression)

class Assign(FunctionAssignment):
    def __str__(self):
        return "%s := %s" % (self.fluent, self.expression)

class Increase(FunctionAssignment):
    pass
