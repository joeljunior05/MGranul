package br.org;

import org.springframework.boot.autoconfigure.*;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.web.bind.annotation.*;

import java.util.HashMap;

@RestController
@EnableAutoConfiguration
public class MGranulServer
{
    @RequestMapping("/")
    public String home() {
        return "Hello World!";
    }

    @GetMapping("/teste")
    public Image greeting() {
        return new Image("Img1", "250 locais");
    }

    public static void main( String[] args )
    {
        HashMap<String, Object> props = new HashMap<String, Object>();

        props.put("server.port", 3418);
        new SpringApplicationBuilder()
                .sources(MGranulServer.class)
                .properties(props)
                .run(args);

    }
}
