from mgkit.io import open_file


def parse_glimmer3(file_handle):
    """
    Parses an ouput file from glimmer3 and yields the header and prediction
    lines. Used to feed the :func:`mgkit.io.gff.from_glimmer3` function.

    Arguments:
        file_handle (str, file): file name or file handle to read from

    Yields:
        tuple: first element is the sequence of the predicted gene and the
        second is the prediction line
    """
    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')

    curr_seq = ''
    predictions = []

    for line in file_handle:
        line = line.strip()
        if line.startswith('>'):
            if len(predictions) > 0:
                for prediction in predictions:
                    yield curr_seq, prediction
            curr_seq = line[1:]
            predictions = []
        else:
            if line != '':
                predictions.append(line)
    else:
        if len(predictions) > 0:
            for prediction in predictions:
                yield curr_seq, prediction
