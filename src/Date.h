#ifndef DATE_H
#define DATE_H

class Date
{
private:
	int m_nMonth;
	int m_nDay;
	int m_nYear;
	
	Date() {}

public:
	Date(int nMonth, int nDay, int nYear);
	void SetDate(int nMonth, int nDay, int nYear);
	int GetMonth() { return m_nMonth;}
	int GetDay() { return m_nDay;}
	int GetYear() {return m_nYear;}
};
#endif
