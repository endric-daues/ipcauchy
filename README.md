![Screenshot](ipcauchy.jpeg)


Welcome! The code contained here follows the theory and numerical study described in (http://www.optimization-online.org/DB_HTML/2020/06/7858.html). Users can download the directory and use the **runner.py** file to count solutions to Knapsack instances with a single equality constraint using a complex path integral that follows a circular, elliptic or shortest path. This code is compatible with Python 3.8.

To install the dependencies use

```
pip install -r requirements.txt
```

To demonstrate the integration methods used in the numerical study, we have provided three example instances in the **instances** directory. You can add your own instances using the simple **json** format. 

```
{
    "name": "instance 1",
    "a": [2, 2, 3, 4, 5],
    "b": 10
}
```

Use the **runner.py** file to run the code with parameters using the following parameters in your terminal.

```
python runner.py --method [integration method] --file [instance json] --N [Optional (SP): Number of Angular Nodes] --r [Optional (SP): radial distance]
```

You can also view the parameter help by using

```
python runner.py --help
```

Here are some examples to get you started.


**Count solutions to Pisinger Instance P1 using a circular integration path**
```
conda run python runner.py --method circle --file instances/pisinger_instance_p1.json
```


**Count solutions to Pisinger Instance P1 using an elliptic integration path**
```
conda run python runner.py --method ellipse --file instances/pisinger_instance_p1.json
```


**Count solutions to Pisinger Instance P1 using a the shortest path integration method (N=36, r=0.001)**
```
conda run python runner.py --method shortest_path --file instances/pisinger_instance_p1.json
```


**Count solutions to Pisinger Instance P20 using a the shortest path integration method (N=360, r=0.001)**
```
conda run python runner.py --method shortest_path --file instances/pisinger_instance_p20.json --N 360 --r 0.001
```