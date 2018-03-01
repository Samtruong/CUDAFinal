#include <cstdlib>
#include <iostream>
using namespace std;

struct node
{
	int info;
	struct node * next;
}*front, *rear;

class queue
{
public:
	void push(int);
	void print();
	struct node* pop();

	queue()
	{
		front = NULL;
		rear = NULL;
	}
};

void queue::push(int info)
{
	node *toInsert;
	toInsert = new struct node;
	toInsert -> info = info;
	toInsert -> next = NULL;
	if (front == NULL)
	{
		front = toInsert;
	}
	else
	{
		rear -> next = toInsert;
	}
	rear = toInsert;
}


struct node * queue::pop()
{
	node *toPop;
	toPop = front;
	front = front -> next;
	return toPop;
}

void queue::print()
{
	node * ptr;
	ptr = front;
	while (ptr != NULL)
	{
		cout << ptr -> info << " ";
		ptr = ptr -> next;
	}
	cout << endl;
}


int main ()
{
	
}