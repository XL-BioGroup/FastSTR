"""Chunked FASTA reader with the installed SimpleFastaParser semantics.

Only column-zero '>' begins a record. TextIO retains the original encoding
and universal-newline behavior. Preserve titles, IDs, empty records, order,
internal tabs and Python slice semantics; no alphabet filtering is added.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_BLOCK_SIZE = 8 * 1024 * 1024
_DELETE_SPACE_NEWLINE = str.maketrans('', '', ' \n')
_OTHER_ASCII_WHITESPACE = '\t\v\f\r\x1c\x1d\x1e\x1f'
_UPPER_ASCII = bytes(range(256)).upper()


def _record_texts(handle, block_size=_BLOCK_SIZE):
    parts = None
    at_line_start = True
    while True:
        chunk = handle.read(block_size)
        if not chunk:
            break
        cursor = 0
        if at_line_start and chunk.startswith('>'):
            position = 0
        else:
            boundary = chunk.find('\n>')
            position = boundary + 1 if boundary >= 0 else -1
        while position >= 0:
            if parts is not None:
                parts.append(chunk[cursor:position])
                record = ''.join(parts)
                parts = []
                yield record
            else:
                parts = []
            cursor = position + 1
            boundary = chunk.find('\n>', cursor)
            position = boundary + 1 if boundary >= 0 else -1
        if parts is not None:
            parts.append(chunk[cursor:])
        at_line_start = chunk.endswith('\n')
    if parts is not None:
        yield ''.join(parts)


def _sequence_text(body):
    # In ordinary ASCII FASTA, only spaces/newlines need removal. Use C-level
    # searches and translate; avoid a Python object for every 80-base line.
    if body.isascii() and not any(c in body for c in _OTHER_ASCII_WHITESPACE):
        return body.translate(_DELETE_SPACE_NEWLINE)
    # SimpleFastaParser rstrips each line, then removes spaces/CR, but keeps
    # INTERNAL tabs. splitlines() or removing all whitespace would differ.
    return ''.join(line.rstrip() for line in body.split('\n')).replace(' ', '').replace('\r', '')


def _sequence_bytes(body):
    if body.isascii() and not any(c in body for c in _OTHER_ASCII_WHITESPACE):
        # Fuse ASCII upper and the exact two-character deletion in one native
        # bytes pass. No alphabet filtering; unknown ASCII symbols survive.
        return body.encode('ASCII').translate(_UPPER_ASCII, b' \n')
    return _sequence_text(body).encode('ASCII').upper()


def read_fasta_reference(path, start_index=1, end_index=0):
    records = []
    with open(path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if end_index == 0:
                record.seq = record.seq.upper()[start_index - 1:]
            else:
                record.seq = record.seq.upper()[start_index - 1:end_index]
            records.append(record)
    return records


def _read_fasta(path, start_index, end_index):
    records = []
    end = None if end_index == 0 else end_index
    with open(path, 'r') as handle:
        for text in _record_texts(handle):
            title, _, body = text.partition('\n')
            title = title.rstrip()
            try:
                first_word = title.split(None, 1)[0]
            except IndexError:
                assert not title, repr(title)
                first_word = ''
            # Validate/normalize the full record before slicing, as before.
            raw = _sequence_bytes(body)
            sequence = Seq(raw[start_index - 1:end])
            records.append(SeqRecord(sequence, id=first_word, name=first_word, description=title))
    return records


def read_fasta(path, start_index=1, end_index=0):
    try:
        return _read_fasta(path, start_index, end_index)
    except UnicodeError:
        # Larger text reads can discover invalid encoding sooner. Replay the
        # legacy parser to retain its failure ordering/messages for such files.
        return read_fasta_reference(path, start_index, end_index)
