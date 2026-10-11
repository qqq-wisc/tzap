"""Numerical conversion of tzap's emitted constant angle expressions."""

from __future__ import annotations

import ast
import math


def angle_radians(source: str) -> float:
    """Convert the restricted constant expression grammar without executing code."""

    def value(node: ast.AST) -> float:
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id == "pi":
            return math.pi
        if isinstance(node, ast.UnaryOp):
            operand = value(node.operand)
            if isinstance(node.op, ast.UAdd):
                return operand
            if isinstance(node.op, ast.USub):
                return -operand
        if isinstance(node, ast.BinOp):
            left, right = value(node.left), value(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError("unsupported angle expression")

    expression = ast.parse(source, mode="eval")
    result = value(expression.body)
    if not math.isfinite(result):
        raise ValueError("angle must be finite")
    return result
