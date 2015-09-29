function formToJSON(formIdentifier) {
    var json = {};

    function addToJSON(json, name, value) {
        if (name.length == 1) {
            json[name[0]] = value;
        }

        else {
            if (json[name[0]] == null) {
                json[name[0]] = {};
            }

            addToJSON(json[name[0]], name.slice(1), value);
        }
    }

    $(document).find("input." + formIdentifier).each(function() {
        if ($(this).attr("name") != "submit") {
            addToJSON(json, $(this).attr("name").split('.'), $(this).val());
        }
    });

    return json;
}

function buildRedirect(data, status, jqXHR) {
    alert("made it");
}

function addFormHandler(url, formIdentifier) {
    $("#submit").click(function() {
        $.ajax({
            type: "POST",
            url: url,
            contentType: "application/json",
            dataType: "json",
            data: JSON.stringify(formToJSON(formIdentifier)),
            success: function(response) {alert(response)}
        });
    });
}

var checkAdminStatus = function() {
    var cookies = document.cookie;
    var lgn = cookies.indexOf("Exologin");

    if (lgn != -1) {
        $("#login").hide();
        $("#logout").show();
    }

    else {
        $("#login").show();
        $("#logout").hide();
    }
};

var main = function() {
    checkAdminStatus();
};

$(document).ready(main());