import argparse
import hashlib
import itertools
import json
from pathlib import Path


def canonical(row):
    def colours(values):
        labels = {}
        return tuple(labels.setdefault(value, len(labels)) for value in values)
    return min(colours(row), colours(tuple(reversed(row))))


def multiply(left, right):
    product = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product[left_index + right_index] += left_value * right_value
    return product


def characteristic(matrix):
    size = len(matrix)
    current = [[int(row == column) for column in range(size)] for row in range(size)]
    coefficients = [1]
    for degree in range(1, size + 1):
        following = [[sum(matrix[row][inner] * current[inner][column]
                          for inner in range(size) if matrix[row][inner])
                      for column in range(size)] for row in range(size)]
        trace = sum(following[index][index] for index in range(size))
        assert trace % degree == 0
        coefficient = -trace // degree
        coefficients.append(coefficient)
        for index in range(size):
            following[index][index] += coefficient
        current = following
    assert not any(value for row in current for value in row)
    while coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def gcd_certificate(left, right):
    for prime in [101, 103, 107, 109, 113]:
        if left[-1] % prime == 0 or right[-1] % prime == 0:
            continue
        before = [value % prime for value in left]
        after = [value % prime for value in right]
        while after:
            while len(before) >= len(after):
                factor = before[-1] * pow(after[-1], -1, prime) % prime
                offset = len(before) - len(after)
                for index, value in enumerate(after):
                    before[index + offset] = (before[index + offset] - factor * value) % prime
                while before and before[-1] == 0:
                    before.pop()
            before, after = after, before
        if len(before) == 1:
            return {'prime': prime, 'last_constant': before[0], 'gcd_degree': 0}
    raise ArithmeticError('No coprimality certificate found')


def construct(size, direction):
    rows = list(itertools.product(range(3), repeat=size))
    classes = sorted({canonical(row) for row in rows})
    positions = {row: index for index, row in enumerate(classes)}
    labelled = [(row, positions[canonical(row)]) for row in rows]
    horizontal = lambda row: sum(row[index] == row[index + 1] for index in range(size - 1)) if direction == 'column' else 0
    zero_classes = [index for index, row in enumerate(classes) if horizontal(row) == 0]
    zero_positions = {index: position for position, index in enumerate(zero_classes)}
    zero_matrix = [[0] * len(zero_classes) for index in zero_classes]
    edges = []
    for before_index, before in enumerate(classes):
        totals = {}
        for after, after_index in labelled:
            cost = horizontal(after) + sum((before[index] == after[index + 1]) +
                    (before[index + 1] == after[index]) for index in range(size - 1))
            if direction == 'row':
                cost += sum(left == right for left, right in zip(before, after))
            if cost <= 2:
                key = after_index, cost
                totals[key] = totals.get(key, 0) + 1
        for (after_index, cost), weight in totals.items():
            if cost == 0 and before_index in zero_positions and after_index in zero_positions:
                zero_matrix[zero_positions[before_index]][zero_positions[after_index]] += weight
            for remaining in range(cost, 3):
                edges.append((3 * before_index + remaining, 3 * after_index + remaining - cost, weight))
    vector = [0] * (3 * len(classes))
    for index, row in enumerate(classes):
        cost = horizontal(row)
        if cost <= 2:
            vector[3 * index + 2 - cost] = sum(after_index == index for after, after_index in labelled)
    return vector, edges, zero_matrix


def verify(size, direction, first, vector_first, expected_order):
    vector, edges, zero_matrix = construct(size, direction)
    base = characteristic(zero_matrix)
    denominator = multiply(multiply(base, base), base)
    order = len(denominator) - 1
    assert order == expected_order
    vectors = [vector]
    for step in range(max(vector_first, 26)):
        following = [0] * len(vector)
        for origin, destination, weight in edges:
            following[destination] += weight * vector[origin]
        vector = following
        vectors.append(vector)
    target = vector_first - 1
    residual = [sum(denominator[lag] * vectors[target - lag][coordinate]
                    for lag in range(order + 1)) for coordinate in range(len(vector))]
    assert not any(residual)
    counts = [0]
    for length, current in enumerate(vectors, 1):
        total = sum(current[::3])
        height, width = (length, size) if direction == 'column' else (size, length)
        constant = (width - 1) * (3 * height - 2) == 2
        assert (total + 3 * constant) % 6 == 0
        counts.append((total + 3 * constant) // 6)
    numerator = [sum(denominator[lag] * counts[index - lag]
                     for lag in range(min(index, order) + 1)) for index in range(len(counts))]
    assert not any(numerator[first:])
    assert numerator[first - 1] != 0
    numerator = numerator[:first]
    certificate = gcd_certificate(numerator, denominator)
    result = {'direction': direction, 'size': size, 'states': len(vector),
            'transitions': len(edges), 'zero_defect_states': len(zero_matrix),
            'zero_defect_denominator_ascending': base,
            'denominator_ascending': denominator, 'numerator_ascending': numerator,
            'minimal_eventual_order': order, 'first_recurrence_index': first,
            'residual_at_previous_index': numerator[-1],
            'full_vector_identity_index': vector_first,
            'identity_coordinates_checked': len(vector), 'nonzero_residual_coordinates': 0,
            'gcd_certificate': certificate, 'initial_counts': counts[1:],
            'vector_table_sha256': hashlib.sha256(json.dumps(vectors, separators=(',', ':')).encode()).hexdigest()}
    if direction == 'column' and size == 5:
        zero_counts = [0] + [sum(current[2::3]) for current in vectors]
        zero_numerator = [sum(base[lag] * zero_counts[index - lag]
                             for lag in range(min(index, len(base) - 1) + 1))
                          for index in range(len(zero_counts))]
        expected = multiply([0, 6], multiply([4, -7, 1], [2, -8, -5, 3]))
        assert zero_numerator[:len(expected)] == expected
        assert not any(zero_numerator[len(expected):])
        assert len(zero_matrix) == 6
        result['zero_defect_companion'] = {'numerator_ascending': expected,
            'denominator_ascending': base, 'initial_labelled_counts': zero_counts[1:],
            'gcd_certificate': gcd_certificate(expected, base)}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--table', type=Path)
    arguments = parser.parse_args()
    specifications = [('column', 2, 6, 6, 3), ('column', 3, 8, 8, 6),
                      ('column', 4, 10, 11, 6), ('column', 5, 21, 21, 18),
                      ('column', 6, 31, 32, 27), ('column', 7, 61, 62, 57),
                      ('row', 2, 7, 7, 6), ('row', 3, 13, 13, 9),
                      ('row', 4, 22, 22, 18), ('row', 5, 34, 34, 30),
                      ('row', 6, 73, 73, 69)]
    results = [verify(size, direction, first, vector_first, order)
               for direction, size, first, vector_first, order in specifications]
    row_two = next(item for item in results if item['direction'] == 'row' and item['size'] == 2)
    assert row_two['denominator_ascending'] == [1, -6, 6, 16, -12, -24, -8]
    columns = {item['size']: item for item in results if item['direction'] == 'column'}
    assert columns[2]['numerator_ascending'] == [0, 0, 6, -15, 13, -2]
    for index, count in enumerate(columns[2]['initial_counts'], 1):
        if index >= 5:
            assert count == 3 ** (index - 5) * (32 * index ** 2 + 31 * index - 30)
    assert columns[2]['initial_counts'][2:4] == [39, 202]
    mathar = [0] + columns[3]['denominator_ascending']
    correction = multiply([0, 0, -4], multiply([-5, 2], [2, -17, 47, -45, 12]))
    mathar += [0] * (len(correction) - len(mathar))
    mathar = [value + correction[index] for index, value in enumerate(mathar)]
    assert mathar == columns[3]['numerator_ascending']
    assert columns[4]['numerator_ascending'] == [0, 3, 59, -976, 4803, -9539, 7326, -1825, -76, 25]
    row_one = [0, 0, 0, 1] + [(width - 1) * (width - 2) // 2 * 2 ** (width - 4) for width in range(4, 31)]
    assert row_one[6] - 6 * row_one[5] + 12 * row_one[4] - 8 * row_one[3] == -4
    assert all(row_one[width] - 6 * row_one[width - 1] + 12 * row_one[width - 2] - 8 * row_one[width - 3] == 0 for width in range(7, 31))
    table_data = {'matched': 0}
    if arguments.table:
        expected = {}
        for line in arguments.table.read_text().splitlines():
            clean = line.split('#')[0].strip()
            if clean:
                index, value = map(int, clean.split())
                expected[index] = value
        sequence_index = 0
        matched = []
        for diagonal in range(2, 30):
            for height in range(1, diagonal):
                width = diagonal - height
                sequence_index += 1
                if sequence_index not in expected:
                    continue
                applicable = [item for item in results
                              if (item['direction'] == 'column' and item['size'] == width and height <= len(item['initial_counts']))
                              or (item['direction'] == 'row' and item['size'] == height and width <= len(item['initial_counts']))]
                observed = []
                for item in applicable:
                    length = height if item['direction'] == 'column' else width
                    observed.append(item['initial_counts'][length - 1])
                if width == 1:
                    observed.append(0)
                if height == 1:
                    observed.append(1 if width == 3 else ((width - 1) * (width - 2) // 2 * 2 ** (width - 4) if width >= 4 else 0))
                if observed:
                    assert all(value == expected[sequence_index] for value in observed), (height, width)
                    matched.append(sequence_index)
        table_data = {'matched': len(matched), 'available': len(expected), 'indices': matched,
                      'sha256': hashlib.sha256(arguments.table.read_bytes()).hexdigest()}
    output = {'schema_version': 1, 'arithmetic': 'exact integers and finite-field Euclidean algorithms',
              'identities': results, 'companion_table': table_data,
              'additional_identities': {'Barker_A229600_generating_function': True,
                  'Barker_A229600_closed_form': True, 'Mathar_A229632_generating_function': True,
                  'Barker_A229633_generating_function': True,
                  'row_one_recurrence_previous_residual': -4}}
    text = json.dumps(output, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
