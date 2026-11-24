import argparse
import hashlib
import itertools
import json
from pathlib import Path

COEFFICIENTS = [8, -228, 2598, -14887, 42804, -40719, -84933, 243513,
                -118011, -230994, 305043, -55215, -116010, 94839, -33687,
                6569, -723, 42]


def multiply(left, right):
    product = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product[left_index + right_index] += left_value * right_value
    return product


def polynomial_gcd_modulo(left, right, prime):
    left = [value % prime for value in left]
    right = [value % prime for value in right]
    degrees = []
    while right:
        while len(left) >= len(right):
            factor = left[-1] * pow(right[-1], -1, prime) % prime
            offset = len(left) - len(right)
            for index, value in enumerate(right):
                left[index + offset] = (left[index + offset] - factor * value) % prime
            while left and left[-1] == 0:
                left.pop()
        degrees.append(len(left) - 1)
        left, right = right, left
    return left, degrees


def sibling_identity(width, coefficients, shift):
    rows = list(itertools.product(range(3), repeat=width))
    states = [(index, remaining) for index in range(len(rows))
              for remaining in range(3)]
    positions = {state: index for index, state in enumerate(states)}
    horizontal = [sum(row[column] == row[column + 1] for column in range(width - 1))
                  for row in rows]
    edges = []
    for before_index, before in enumerate(rows):
        for after_index, after in enumerate(rows):
            cost = horizontal[after_index] + sum(
                (before[column] == after[column + 1]) +
                (before[column + 1] == after[column]) for column in range(width - 1))
            for remaining in range(cost, 3):
                edges.append((positions[before_index, remaining],
                              positions[after_index, remaining - cost]))
    vector = [0] * len(states)
    for row_index, cost in enumerate(horizontal):
        if cost <= 2:
            vector[positions[row_index, 2 - cost]] = 1
    vectors = [vector]
    for step in range(shift + len(coefficients)):
        following = [0] * len(states)
        for origin, destination in edges:
            following[destination] += vector[origin]
        vector = following
        vectors.append(vector)
    residual = [vectors[-1][coordinate] - sum(
        coefficient * vectors[shift + index][coordinate]
        for index, coefficient in enumerate(coefficients))
        for coordinate in range(len(states))]
    if any(residual):
        raise ArithmeticError(f'The width-{width} transfer-vector identity is false')
    return {'width': width, 'coefficients_ascending': coefficients,
            'vector_shift': shift, 'identity_coordinates_checked': len(states),
            'nonzero_residual_coordinates': 0,
            'recurrence_first_index': shift + len(coefficients) + 1}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--table', type=Path)
    args = parser.parse_args()
    rows = list(itertools.product(range(3), repeat=5))
    horizontal = [sum(row[column] == row[column + 1] for column in range(4))
                  for row in rows]
    states = [(index, remaining) for index in range(len(rows))
              for remaining in range(3)]
    positions = {state: index for index, state in enumerate(states)}
    edges = []
    for before_index, before in enumerate(rows):
        for after_index, after in enumerate(rows):
            cost = horizontal[after_index] + sum(
                (before[column] == after[column + 1]) +
                (before[column + 1] == after[column]) for column in range(4))
            for remaining in range(cost, 3):
                edges.append((positions[before_index, remaining],
                              positions[after_index, remaining - cost]))
    vector = [0] * len(states)
    for row_index, cost in enumerate(horizontal):
        if cost <= 2:
            vector[positions[row_index, 2 - cost]] = 1
    vectors = [vector]
    for step in range(20):
        following = [0] * len(states)
        for origin, destination in edges:
            following[destination] += vector[origin]
        vector = following
        vectors.append(vector)
    residual = [vectors[20][coordinate] - sum(
        coefficient * vectors[2 + shift][coordinate]
        for shift, coefficient in enumerate(COEFFICIENTS))
        for coordinate in range(len(states))]
    if any(residual):
        raise ArithmeticError('The transfer-vector identity is false')
    counts = []
    for current in vectors:
        labelled = sum(current[positions[row_index, 0]]
                       for row_index in range(len(rows)))
        if labelled % 6:
            raise ArithmeticError('Colour-orbit divisibility failed')
        counts.append(labelled // 6)
    data = {'schema_version': 1, 'arithmetic': 'exact integers',
            'rows': len(rows), 'states': len(states), 'transitions': len(edges),
            'identity': 'w_20 = sum(c_j * w_(2+j), j=0..17)',
            'coefficients_ascending': COEFFICIENTS,
            'identity_coordinates_checked': len(residual),
            'nonzero_residual_coordinates': sum(value != 0 for value in residual),
            'initial_counts_n_1_to_21': counts,
            'vector_table_sha256': hashlib.sha256(json.dumps(
                vectors, separators=(',', ':')).encode()).hexdigest()}
    denominator = [1] + [-value for value in reversed(COEFFICIENTS)]
    series = [0] + counts
    numerator = [sum(denominator[shift] * series[index - shift]
                     for shift in range(min(index, 18) + 1)) for index in range(21)]
    base = [1, -14, 45, -15, -36, 19, -2]
    if multiply(multiply(base, base), base) != denominator:
        raise ArithmeticError('The denominator factorization is false')
    divisor, degrees = polynomial_gcd_modulo(numerator, denominator, 101)
    if len(divisor) != 1 or numerator[-1] % 101 == 0 or denominator[-1] % 101 == 0:
        raise ArithmeticError('The reduced-denominator certificate failed')
    if numerator[-1] != -80:
        raise ArithmeticError('The recurrence boundary certificate failed')
    data['generating_function'] = {
        'numerator_ascending': numerator, 'denominator_ascending': denominator,
        'denominator_cube_root_ascending': base, 'gcd_prime': 101,
        'gcd_remainder_degrees': degrees, 'gcd_last_nonzero_constant': divisor[0],
        'minimal_eventual_recurrence_order': 18, 'recurrence_residual_at_n_20': -80}
    data['sibling_identities'] = [sibling_identity(2, [27, -27, 9], 2),
                                  sibling_identity(3, [-8, 60, -162, 185, -81, 15], 1)]
    if args.table:
        expected = {}
        for line in args.table.read_text().splitlines():
            clean = line.strip()
            if clean and not clean.startswith('#'):
                index, value = map(int, clean.split())
                expected[index] = value
        matched = 0
        current = vectors[-1]
        for index in range(1, max(expected) + 1):
            if index <= len(counts):
                observed = counts[index - 1]
            else:
                following = [0] * len(states)
                for origin, destination in edges:
                    following[destination] += current[origin]
                current = following
                labelled = sum(current[positions[row_index, 0]]
                               for row_index in range(len(rows)))
                if labelled % 6:
                    raise ArithmeticError('Colour-orbit divisibility failed')
                observed = labelled // 6
            if index in expected:
                if observed != expected[index]:
                    raise ArithmeticError(f'Table mismatch at n={index}')
                matched += 1
        data['tabulated_counts_matched'] = matched
        data['table_sha256'] = hashlib.sha256(args.table.read_bytes()).hexdigest()
    text = json.dumps(data, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
