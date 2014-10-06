/*
 * Created by: Leslie Sanford
 * 
 * Contact: jabberdabber@hotmail.com
 * 
 * Last modified: 09/26/2005
 */

using System;
using System.Collections;

namespace CNNWB.Common
{
	/// <summary>
	/// Represents a simple double-ended-queue collection of objects.
	/// </summary>
	[Serializable()]
	public class Deque : ICollection, IEnumerable, ICloneable
	{
        #region Deque Members

        #region Fields

        // The node at the front of the deque.
        private Node front = null;

        // The node at the back of the deque.
        private Node back = null;

        // The number of elements in the deque.
        private int count = 0;

        // The version of the deque.
        private long version = 0;

        #endregion

        #region Construction

        /// <summary>
        /// Initializes a new instance of the Deque class.
        /// </summary>
		public Deque()
		{
        }

        /// <summary>
        /// Initializes a new instance of the Deque class that contains 
        /// elements copied from the specified collection.
        /// </summary>
        /// <param name="col">
        /// The ICollection to copy elements from.
        /// </param>
        public Deque(ICollection col)
        {
            #region Preconditions

            if(col == null)
            {
                throw new ArgumentNullException("col");
            }

            #endregion

            foreach(object obj in col)
            {
                PushBack(obj);
            }
        }

        #endregion

        #region Methods

        /// <summary>
        /// Removes all objects from the Deque.
        /// </summary>
        public virtual void Clear()
        {
            count = 0;

            front = back = null;

            version++;
        }

        /// <summary>
        /// Determines whether or not an element is in the Deque.
        /// </summary>
        /// <param name="obj">
        /// The Object to locate in the Deque.
        /// </param>
        /// <returns>
        /// <b>true</b> if <i>obj</i> if found in the Deque; otherwise, 
        /// <b>false</b>.
        /// </returns>
        public virtual bool Contains(object obj)
        {
            foreach(object o in this)
            {
                if(o == null && obj == null)
                {
                    return true;
                }
                else if(o.Equals(obj))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Inserts an object at the front of the Deque.
        /// </summary>
        /// <param name="obj">
        /// The object to push onto the deque;
        /// </param>
        public virtual void PushFront(object obj)
        {
            // The new node to add to the front of the deque.
            Node n = new Node(obj);

            // Link the new node to the front node. The current front node at 
            // the front of the deque is now the second node in the deque.
            n.Next = front;

            // If the deque isn't empty
            if(Count > 0)
            {
                // Link the current front to the new node.
                front.Previous = n;
            }

            // Make the new node the front of the deque.
            front = n;            

            // Keep track of the number of elements in the deque.
            count++;

            // If this is the first element in the deque.
            if(Count == 1)
            {
                // The front and back nodes are the same.
                back = front;
            }

            version++;
        }

        /// <summary>
        /// Inserts an object at the back of the Deque.
        /// </summary>
        /// <param name="obj">
        /// The object to push onto the deque;
        /// </param>
        public virtual void PushBack(object obj)
        {
            // The new node to add to the back of the deque.
            Node n = new Node(obj);
            
            // Link the new node to the back node. The current back node at 
            // the back of the deque is now the second to the last node in the
            // deque.
            n.Previous = back;

            // If the deque is not empty.
            if(Count > 0)
            {
                // Link the current back node to the new node.
                back.Next = n;
            }

            // Make the new node the back of the deque.
            back = n;            

            // Keep track of the number of elements in the deque.
            count++;

            // If this is the first element in the deque.
            if(Count == 1)
            {
                // The front and back nodes are the same.
                front = back;
            }

            version++;
        }        

        /// <summary>
        /// Removes and returns the object at the front of the Deque.
        /// </summary>
        /// <returns>
        /// The object at the front of the Deque.
        /// </returns>
        /// <exception cref="InvalidOperationException">
        /// The Deque is empty.
        /// </exception>
        public virtual object PopFront()
        {
            #region Preconditions

            if(Count == 0)
            {
                throw new InvalidOperationException("Deque is empty.");
            }

            #endregion

            // Get the object at the front of the deque.
            object obj = front.Value;

            // Move the front back one node.
            front = front.Next;

            // Keep track of the number of nodes in the deque.
            count--;

            // If the deque is not empty.
            if(Count > 0)
            {
                // Tie off the previous link in the front node.
                front.Previous = null;
            }
            // Else the deque is empty.
            else
            {
                // Indicate that there is no back node.
                back = null;
            }           

            version++;

            return obj;
        }

        /// <summary>
        /// Removes and returns the object at the back of the Deque.
        /// </summary>
        /// <returns>
        /// The object at the back of the Deque.
        /// </returns>
        /// <exception cref="InvalidOperationException">
        /// The Deque is empty.
        /// </exception>
        public virtual object PopBack()
        {
            #region Preconditions

            if(Count == 0)
            {
                throw new InvalidOperationException("Deque is empty.");
            }

            #endregion

            // Get the object at the back of the deque.
            object obj = back.Value;

            // Move back node forward one node.
            back = back.Previous;

            // Keep track of the number of nodes in the deque.
            count--;

            // If the deque is not empty.
            if(Count > 0)
            {
                // Tie off the next link in the back node.
                back.Next = null;
            }
            // Else the deque is empty.
            else
            {
                // Indicate that there is no front node.
                front = null;
            }

            version++;

            return obj;
        }

        /// <summary>
        /// Returns the object at the front of the Deque without removing it.
        /// </summary>
        /// <returns>
        /// The object at the front of the Deque.
        /// </returns>
        /// <exception cref="InvalidOperationException">
        /// The Deque is empty.
        /// </exception>
        public virtual object PeekFront()
        {
            #region Preconditions

            if(Count == 0)
            {
                throw new InvalidOperationException("Deque is empty.");
            }

            #endregion

            return front.Value;
        }

        /// <summary>
        /// Returns the object at the back of the Deque without removing it.
        /// </summary>
        /// <returns>
        /// The object at the back of the Deque.
        /// </returns>
        /// <exception cref="InvalidOperationException">
        /// The Deque is empty.
        /// </exception>
        public virtual object PeekBack()
        {
            #region Preconditions

            if(Count == 0)
            {
                throw new InvalidOperationException("Deque is empty.");
            }

            #endregion

            return back.Value;
        }

        /// <summary>
        /// Copies the Deque to a new array.
        /// </summary>
        /// <returns>
        /// A new array containing copies of the elements of the Deque.
        /// </returns>
        public virtual object[] ToArray()
        {
            object[] array = new object[Count];
            int index = 0;

            foreach(object obj in this)
            {
                array[index] = obj;
                index++;
            }

            return array;
        }

        /// <summary>
        /// Returns a synchronized (thread-safe) wrapper for the Deque.
        /// </summary>
        /// <param name="deque">
        /// The Deque to synchronize.
        /// </param>
        /// <returns>
        /// A synchronized wrapper around the Deque.
        /// </returns>
        public static Deque Synchronized(Deque deque)
        {            
            #region Preconditions

            if(deque == null)
            {
                throw new ArgumentNullException("deque");
            }

            #endregion

            return new SynchronizedDeque(deque);
        }        

        #endregion

        #region Node Class

        // Represents a node in the deque.
        [Serializable()]
        private class Node
        {
            private object value;

            private Node previous = null;

            private Node next = null;

            public Node(object value)
            {
                this.value = value;
            }

            public object Value
            {
                get
                {
                    return value;
                }
            }

            public Node Previous
            {
                get
                {
                    return previous;
                }
                set
                {
                    previous = value;
                }
            }

            public Node Next
            {
                get
                {
                    return next;
                }
                set
                {
                    next = value;
                }
            }
        }

        #endregion

        #region DequeEnumerator Class

        [Serializable()]
        private class DequeEnumerator : IEnumerator
        {
            private Deque owner;

            private Node current;

            private bool beforeBeginning = true;

            private long version;

            public DequeEnumerator(Deque owner)
            {
                this.owner = owner;
                current = owner.front;
                this.version = owner.version;
            }

            #region IEnumerator Members

            public void Reset()
            {
                #region Preconditions

                if(version != owner.version)
                {
                    throw new InvalidOperationException(
                        "The Deque was modified after the enumerator was created.");
                }

                #endregion

                beforeBeginning = true;
            }

            public object Current
            {
                get
                {
                    #region Preconditions

                    if(beforeBeginning || current == null)
                    {
                        throw new InvalidOperationException(
                            "The enumerator is positioned before the first " +
                            "element of the Deque or after the last element.");
                    }

                    #endregion

                    return current.Value;
                }
            }

            public bool MoveNext()
            {
                #region Preconditions

                if(version != owner.version)
                {
                    throw new InvalidOperationException(
                        "The Deque was modified after the enumerator was created.");
                }

                #endregion

                bool result = false;

                // If the enumerator is positioned before the front of the 
                // deque.
                if(beforeBeginning)
                {
                    // Position the enumerator at the front of the deque.
                    current = owner.front;

                    // Indicate that the enumerator is no longer positioned
                    // before the beginning of the deque.
                    beforeBeginning = false;
                }
                // Else if the enumerator has not yet reached the end of the 
                // deque.
                else if(current != null)
                {
                    // Move to the next element in the deque.
                    current = current.Next;
                }

                // If the enumerator has not reached the end of the deque.
                if(current != null)
                {
                    // Indicate that the enumerator has not reached the end of
                    // the deque.
                    result = true;
                }

                return result;
            }

            #endregion
        }

        #endregion

        #region SynchronizedDeque Class

        // Implements a synchronization wrapper around a deque.
        [Serializable()]
        private class SynchronizedDeque : Deque
        {
            #region SynchronziedDeque Members

            #region Fields

            // The wrapped deque.
            private Deque deque;

            // The object to lock on.
            private object root;

            #endregion

            #region Construction

            public SynchronizedDeque(Deque deque)
            {
                this.deque = deque;
                this.root = deque.SyncRoot;
            }

            #endregion

            #region Methods

            public override void Clear()
            {
                lock(root)
                {
                    deque.Clear();
                }
            }

            public override bool Contains(object obj)
            {
                bool result;

                lock(root)
                {
                    result = deque.Contains(obj);
                }

                return result;
            }

            public override void PushFront(object obj)
            {
                lock(root)
                {
                    deque.PushFront(obj);
                }
            }

            public override void PushBack(object obj)
            {
                lock(root)
                {
                    deque.PushBack(obj);
                }
            }

            public override object PopFront()
            {
                object obj;

                lock(root)
                {
                    obj = deque.PopFront();
                }
            
                return obj;
            }

            public override object PopBack()
            {
                object obj;

                lock(root)
                {
                    obj = deque.PopBack();
                }

                return obj;
            }

            public override object PeekFront()
            {
                object obj;

                lock(root)
                {
                    obj = deque.PeekFront();
                }

                return obj;
            }

            public override object PeekBack()
            {
                object obj;

                lock(root)
                {
                    obj = deque.PeekBack();
                }

                return obj;
            }

            public override object[] ToArray()
            {
                object[] array;

                lock(root)
                {
                    array = deque.ToArray();
                }

                return array;
            }

            public override object Clone()
            {
                object clone;

                lock(root)
                {
                    clone = deque.Clone();
                }

                return clone;
            }

            public override void CopyTo(Array array, int index)
            {
                lock(root)
                {
                    deque.CopyTo(array, index);
                }
            }

            public override IEnumerator GetEnumerator()
            {
                IEnumerator e;

                lock(root)
                {
                    e = deque.GetEnumerator();
                }

                return e;
            }

            #endregion

            #region Properties

            public override int Count
            {
                get
                {
                    int count;

                    lock(root)
                    {
                        count = deque.Count;
                    }

                    return count;
                }
            }

            public override bool IsSynchronized
            {
                get
                {
                    return true;
                }
            }

            #endregion

            #endregion
        }

        #endregion

        #endregion

        #region ICollection Members

        /// <summary>
        /// Gets a value indicating whether access to the Deque is synchronized 
        /// (thread-safe).
        /// </summary>
        public virtual bool IsSynchronized
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets the number of elements contained in the Deque.
        /// </summary>
        public virtual int Count
        {
            get
            {
                return count;
            }
        }

        /// <summary>
        /// Copies the Deque elements to an existing one-dimensional Array, 
        /// starting at the specified array index.
        /// </summary>
        /// <param name="array">
        /// The one-dimensional Array that is the destination of the elements 
        /// copied from Deque. The Array must have zero-based indexing. 
        /// </param>
        /// <param name="index">
        /// The zero-based index in array at which copying begins. 
        /// </param>
        public virtual void CopyTo(Array array, int index)
        {
            #region Preconditions

            if(array == null)
            {
                throw new ArgumentNullException("array");
            }
            else if(index < 0)
            {
                throw new ArgumentOutOfRangeException("index", index,
                    "Index is less than zero.");
            }
            else if(array.Rank > 1)
            {
                throw new ArgumentException("Array is multidimensional.");
            }
            else if(index >= array.Length)
            {
                throw new ArgumentException("Index is equal to or greater " +
                    "than the length of array.");
            }
            else if(Count > array.Length - index)
            {
                throw new ArgumentException(
                    "The number of elements in the source Deque is greater " +
                    "than the available space from index to the end of the " +
                    "destination array.");
            }

            #endregion

            int i = index;

            foreach(object obj in this)
            {
                array.SetValue(obj, i);
                i++;
            }
        }

        /// <summary>
        /// Gets an object that can be used to synchronize access to the Deque.
        /// </summary>
        public virtual object SyncRoot
        {
            get
            {
                return this;
            }
        }
        
        #endregion

        #region IEnumerable Members

        /// <summary>
        /// Returns an enumerator that can iterate through the Deque.
        /// </summary>
        /// <returns>
        /// An IEnumerator for the Deque.
        /// </returns>
        public virtual IEnumerator GetEnumerator()
        {
            return new DequeEnumerator(this);
        }

        #endregion

        #region ICloneable Members

        /// <summary>
        /// Creates a shallow copy of the Deque.
        /// </summary>
        /// <returns>
        /// A shalloe copy of the Deque.
        /// </returns>
        public virtual object Clone()
        {
            Deque clone = new Deque(this);

            clone.version = this.version;

            return clone;
        }

        #endregion
    }
}
