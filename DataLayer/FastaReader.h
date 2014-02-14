#ifndef FASTAREADER_H
#define FASTAREADER_H 1

#include "Sequence.h"
#include "StringUtil.h" // for chomp
#include <cassert>
#include <cstdlib> // for exit
#include <fstream>
#include <istream>
#include <iostream> // for debugging, remove before committing
#include <limits> // for numeric_limits
#include <ostream>
#include <algorithm>

#include <cstdio>
#include <cstring>

static inline int fpeek(FILE * stream)
{
	int c;
	c = fgetc(stream);
	ungetc(c, stream);
	return c;
}

/** Read a FASTA, FASTQ, export, qseq or SAM file. */
class FastaReader {
public:
	enum {
		/** Fold lower-case characters to upper-case. */
		FOLD_CASE = 0, NO_FOLD_CASE = 1,
		/** Convert to standard quality. */
		NO_CONVERT_QUALITY = 0, CONVERT_QUALITY = 2,
	};
	bool flagFoldCase()
	{
		return ~m_flags & NO_FOLD_CASE;
	}
	bool flagConvertQual()
	{
		return m_flags & CONVERT_QUALITY;
	}

	FastaReader(const char* path, int flags, int len = 0);

	~FastaReader()
	{
		if (!eof()) {
			std::string line;
			getline(line);
			die() << "expected end-of-file near\n" << line << '\n';
			exit(EXIT_FAILURE);
		}
		fclose(m_in);
		delete m_buff;
	}

	Sequence read(std::string& id, std::string& comment, char& anchor,
			std::string& qual);

	/** Return whether this stream is at end-of-file. */
	bool eof() const
	{
		return m_bstart >= m_bend && feof(m_in);
	}
	;

	/** Return whether this stream is good. */
	operator bool() const
	{
		return !m_fail;
	}

	/** Return the next character of this stream. */
	int peek()
	{
		if (fill_buff())
			return m_buff[m_bstart];
		return EOF;
	}

	/** Interface for manipulators. */
	//FastaReader& operator>>(std::istream& (*f)(std::istream&))
	//{
	//	f(m_in);
	//	return *this;
	//}
	/** Returns the number of unchaste reads. */
	unsigned unchaste() const
	{
		return m_unchaste;
	}

	FastaReader& operator >>(Sequence& seq)
	{
		std::string id, comment, qual;
		char anchor;
		seq = this->read(id, comment, anchor, qual);
		return *this;
	}

private:
	bool fill_buff()
	{
		if (m_bstart >= m_bend) {
			m_bstart = 0;
			m_bend = fread(m_buff, sizeof(char), m_blen, m_in);
		}
		return m_bend != 0;
	}

	/** Read a single line. */
	bool getline(std::string& s)
	{
		s.clear();
		while (!eof()) {
			if (!fill_buff())
				break;
			// find first new line char
			char * end = strchr(m_buff + m_bstart, '\n');
			if (end == NULL || m_buff + m_bend < end)
				end = m_buff + m_bend;
			assert(end >= m_buff);
			s.append(m_buff + m_bstart, end);
			m_bstart = end - m_buff + 1;
			if (m_bstart <= m_bend)
				break;
		}
		chomp(s, '\r');
		m_line++;
		if (s.empty())
			m_fail = true;
		return !m_fail;
	}

	void clear()
	{
		m_fail = false;
	}

	/** Ignore the specified number of lines. */
	bool ignoreLines(unsigned n)
	{
		std::string s;
		bool good = true;
		for (unsigned i = 0; i < n; ++i) {
			if ((good = getline(s)))
				m_line++;
		}
		return good;
	}

	std::ostream& die();
	bool isChaste(const std::string& s, const std::string& line);
	void checkSeqQual(const std::string& s, const std::string& q);

	const char* m_path;
	//std::ifstream m_fin;
	//std::istream& m_in;
	FILE * m_in;

	size_t m_blen, m_bstart, m_bend;
	char * m_buff;
	bool m_fail;

	/** Flags indicating parsing options. */
	int m_flags;

	/** Number of lines read. */
	unsigned m_line;

	/** Count of unchaste reads. */
	unsigned m_unchaste;

	/** Position of the end of the current section. */
	std::streampos m_end;

	/** Trim sequences to this length. 0 is unlimited. */
	const int m_maxLength;
};

/** A FASTA record. */
struct FastaRecord {
	/** Identifier */
	std::string id;
	/** Comment following the first white-space of the header */
	std::string comment;
	/** Anchor base for a colour-space sequence */
	char anchor;
	/** The sequence */
	Sequence seq;

	FastaRecord()
	{
	}
	FastaRecord(const std::string& id, const std::string& comment,
			const Sequence& seq) :
			id(id), comment(comment), anchor(0), seq(seq)
	{
	}

	operator Sequence() const
	{
		return seq;
	}

	FastaRecord& operator=(const std::string& s)
	{
		seq = s;
		return *this;
	}

	size_t size() const
	{
		return seq.size();
	}

	friend FastaReader& operator >>(FastaReader& in, FastaRecord& o)
	{
		std::string q;
		o.seq = in.read(o.id, o.comment, o.anchor, q);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out, const FastaRecord& o)
	{
		out << '>' << o.id;
		if (!o.comment.empty())
			out << ' ' << o.comment;
		return out << '\n' << o.seq << '\n';
	}
};

/** A FASTQ record. */
struct FastqRecord: FastaRecord {
	/** Quality */
	std::string qual;

	FastqRecord()
	{
	}
	FastqRecord(const std::string& id, const std::string& comment,
			const Sequence& seq, const std::string& qual) :
			FastaRecord(id, comment, seq), qual(qual)
	{
		assert(seq.length() == qual.length());
	}

	friend FastaReader& operator >>(FastaReader& in, FastqRecord& o)
	{
		o.seq = in.read(o.id, o.comment, o.anchor, o.qual);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out, const FastqRecord& o)
	{
		if (o.qual.empty())
			return out << static_cast<const FastaRecord&>(o);
		out << '@' << o.id;
		if (!o.comment.empty())
			out << ' ' << o.comment;
		return out << '\n' << o.seq << "\n"
				"+\n" << o.qual << '\n';
	}
};

#endif //FASTAREADER_H
