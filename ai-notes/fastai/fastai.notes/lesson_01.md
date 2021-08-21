
### Lesson 01
* this course starts with coding and doing (learning how) and works towards theory (learning why).

* still a lot of artisanship in ML

* we need to have a good intuitive understanding of what the data looks like and what the labels are, before we can attempt to solve any problems or make any predictions

* some of the first data we'll download is the pets dataset (cited in the nb); 37 breeds of dogs and cats; they used to use a dataset that was made to build classifiers that can distinguish between dog *or* cat. that was a hard problem at the time. now it is much too easy. this is the **fine-grain classification** problem to the easier **coarse-grained classification**.

* how do we get the labels? in ML, **labels** refer to the thing that we are trying to predict. in this case, the labels are in the file names. i.e.,
`.fastai/data/oxford-iiit-pet/images/Siamese_169.jpg`

`ImageDataBunch / DataBunchObject` contains all elements that are needed to build a model: training, validation, and test sets; labels, tabular data, etc.

* we format our image data in a square - we need to pass the same regularized format to the GPU each time to enable the fast computing.
* `size = 224` generally works for most things.

* By default, fastai it will do center-cropping of the input images in the dataset. This combination of cropping and re-sizing enables the normalization of datatypes across all images. This is done slightly differently on purpose to encourage randomization / eliminate systematic error.

* we're already ready to train a model and we'll start with a convolutional neural net (resnet34). This will use a single hidden layer as a classifier. more importantly for now, it will have 37 outputs in which it predicts the probability that the input is one of the 37 potential breed labels.

* resnet34, when it starts, downloads a pre-trained model that was trained on image net - it already knew how to classify general images pretty well (ImageNet) - we're taking something that already does something pretty well and teaching it to do *our thing* really well. When we do this, in general it takes 1/100th of the time to train models with < 1/100th of the data. This is really important!

* we need to avoid overfitting - to do this, we use a validation set that our model never got to see. When we created the `dataBunch`, fastai did this for us automatically. if you're not using a validation set, you don't know if you're overfitting.

* it's normal to have some categories that the model doesn't do very well predicting; these would be revealed in a confusion matrix or in the figure where we show the most inaccurately predicted labels.

* in learning this material, it is important to pay attention to the code notebooks and to what is going into and coming out of the models, rather than looking up theory, etc. the notebooks are key!

* we were able to do this really quickly because of fastai. use the fastai documentation to learn the fastai software deeply. any time they can pick a default, they pick it for you, making it less coding and easier to try things.

* Google Brain and OpenAI are two of the top places to do AI research.

* edge computing - run something on a mobile phone or small device.

* ResNet50 is common and works well across a wide range of applications.

* in fastai, if you have many classes, a confusion matrix can be difficult to read - fastai has a .most_confused() function that prints the most confused combinations of categories by your model.

* each layer of a convnet might find various aspects that differentiate pictures (gradients, shapes, edges, etc.) - eventually this adds up to larger structures (combining previous layers), such that you can compare more complex ideas like "fluffy" vs "smooth" so that by the time you get to the final layer, you can classify your real categories (one breed vs another).
